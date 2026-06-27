"""Sequence diversity analysis for AMP candidate panels.

Two concerns before synthesis:
1. Redundancy — candidates >60% similar to another in the panel occupy the same
   structural space; testing both adds cost without proportional information gain.
2. Family-level risk — if all variants of a scaffold share a structural liability
   (e.g., 5 trypsin sites in all SEED-005 14-mers), the whole family may fail
   together, wiping out the corresponding synthesis investment.

Key outputs:
- Pairwise similarity matrix (Levenshtein-based, no alignment required)
- Cluster assignments (single-linkage at configurable threshold)
- Recommended minimal diverse panel: best 1 candidate per cluster
- Family-level structural warnings
"""
from __future__ import annotations


def levenshtein_similarity(a: str, b: str) -> float:
    """Sequence similarity in [0, 1] via normalised Levenshtein distance.

    Defined as 1 - edit_distance / max(len(a), len(b)).
    Pure Python, O(|a| * |b|) — fine for short peptides.
    """
    a, b = a.upper(), b.upper()
    la, lb = len(a), len(b)
    if la == 0 and lb == 0:
        return 1.0
    dp = list(range(lb + 1))
    for ca in a:
        ndp = [dp[0] + 1]
        for j, cb in enumerate(b):
            ndp.append(min(dp[j + 1] + 1, ndp[j] + 1, dp[j] + (0 if ca == cb else 1)))
        dp = ndp
    return 1.0 - dp[lb] / max(la, lb)


def pairwise_similarity_matrix(
    sequences: list[str],
) -> list[list[float]]:
    """Return an n×n symmetric similarity matrix."""
    n = len(sequences)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        mat[i][i] = 1.0
        for j in range(i + 1, n):
            sim = levenshtein_similarity(sequences[i], sequences[j])
            mat[i][j] = sim
            mat[j][i] = sim
    return mat


def cluster_panel(
    candidates: list[dict],
    similarity_threshold: float = 0.60,
) -> list[dict]:
    """Assign cluster IDs using greedy single-linkage clustering.

    Each candidate dict must have 'candidate_id' and 'sequence'.
    Returns candidates with an added 'cluster_id' field.
    Earlier (higher-priority) candidates seed clusters.
    """
    sequences = [c["sequence"] for c in candidates]
    cluster_ids: list[int] = [-1] * len(candidates)
    next_cluster = 0

    for i, seq_i in enumerate(sequences):
        if cluster_ids[i] != -1:
            continue
        cluster_ids[i] = next_cluster
        for j in range(i + 1, len(sequences)):
            if cluster_ids[j] != -1:
                continue
            if levenshtein_similarity(seq_i, sequences[j]) >= similarity_threshold:
                cluster_ids[j] = next_cluster
        next_cluster += 1

    return [
        {**c, "cluster_id": cluster_ids[i]}
        for i, c in enumerate(candidates)
    ]


def recommend_minimal_diverse_panel(
    clustered: list[dict],
    n_per_cluster: int = 1,
) -> list[dict]:
    """Pick the top `n_per_cluster` candidates per cluster by list order.

    Assumes the input is already sorted by priority (highest first).
    Returns the minimal set that covers all structural families.
    """
    seen: dict[int, int] = {}
    minimal: list[dict] = []
    for c in clustered:
        cid = c["cluster_id"]
        count = seen.get(cid, 0)
        if count < n_per_cluster:
            minimal.append(c)
            seen[cid] = count + 1
    return minimal


def diversity_stats(clustered: list[dict]) -> dict:
    """Compute diversity statistics for a clustered panel."""
    from collections import Counter
    cluster_sizes = Counter(c["cluster_id"] for c in clustered)
    n_clusters = len(cluster_sizes)
    n_candidates = len(clustered)
    redundant = sum(1 for c in clustered if cluster_sizes[c["cluster_id"]] > 1)
    singletons = sum(1 for s in cluster_sizes.values() if s == 1)
    largest_cluster = max(cluster_sizes.values())
    sequences = [c["sequence"] for c in clustered]
    mat = pairwise_similarity_matrix(sequences)
    # Mean off-diagonal similarity
    n = len(sequences)
    if n > 1:
        total = sum(mat[i][j] for i in range(n) for j in range(n) if i != j)
        mean_sim = total / (n * (n - 1))
    else:
        mean_sim = 1.0
    return {
        "n_candidates": n_candidates,
        "n_clusters": n_clusters,
        "n_redundant": redundant,
        "n_singletons": singletons,
        "largest_cluster_size": largest_cluster,
        "mean_pairwise_similarity": round(mean_sim, 4),
        "diversity_score": round(1.0 - mean_sim, 4),
    }


def family_structural_warnings(
    candidates: list[dict],
    min_family_size: int = 3,
) -> list[dict]:
    """Detect families where ALL members share a structural liability.

    Expects each candidate dict to optionally include:
    - 'seed': e.g. 'SEED-005'
    - 'trypsin_sites' (list) from presynth QC
    - 'mu_h' (float)
    - 'methionine_count' (int)

    Returns a list of warning dicts.
    """
    from collections import defaultdict
    families: dict[str, list[dict]] = defaultdict(list)
    for c in candidates:
        seed = c.get("seed", "UNKNOWN")
        families[seed].append(c)

    warnings = []
    for seed, members in families.items():
        if len(members) < min_family_size:
            continue
        # Trypsin sites warning: all members have >= 4 sites
        trypsin_counts = [len(c.get("trypsin_sites", [])) for c in members if "trypsin_sites" in c]
        if trypsin_counts and all(t >= 4 for t in trypsin_counts):
            warnings.append({
                "family": seed,
                "n_members": len(members),
                "warning_type": "TRYPSIN_STABILITY",
                "message": (
                    f"All {len(members)} {seed} members have ≥4 trypsin sites — "
                    "serum stability < 2h expected for the ENTIRE family. "
                    "If goal is in-vivo use, consider N-methylation or D-amino acid substitutions."
                ),
            })
        # High μH warning: all members μH > 0.65
        mu_h_vals = [c.get("mu_h", 0.0) for c in members if "mu_h" in c]
        if mu_h_vals and all(m > 0.65 for m in mu_h_vals):
            warnings.append({
                "family": seed,
                "n_members": len(members),
                "warning_type": "HEMOLYTIC_FAMILY_RISK",
                "message": (
                    f"All {len(members)} {seed} members have μH > 0.65 — "
                    "entire family is in the moderate-to-high hemolytic risk tier. "
                    "Test hemolysis before MIC; if HC50 < 2×MIC, family may lack a therapeutic window."
                ),
            })
    return warnings
