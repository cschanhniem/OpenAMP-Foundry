"""Retrospective AUROC benchmark against known AMPs vs composition-matched decoys.

Tests whether our scoring model ranks known antimicrobial peptides higher than
randomly shuffled decoys with the same amino acid composition.

Benchmark design:
  Positives — 44 sequences from amp_curated_references.csv (confirmed literature AMPs)
  Negatives — 44 per-sequence shuffled decoys (RNG seed=42, same composition)

Key metric: AUROC (area under receiver operating characteristic curve)
  = P(random AMP scores higher than random decoy)
  < 0.55 → model is near-random; do NOT proceed to synthesis
  0.55–0.70 → model has weak but real signal; treat $10k budget with caution
  > 0.70 → model has meaningful discriminative power; proceed to synthesis

IMPORTANT: This benchmark ONLY tests discrimination of order-dependent features
(hydrophobic moment) since composition-based features are identical by design.
A high AUROC here means the amphipathicity signal is real; it does NOT test
whether our nominees are actually antimicrobial.
"""
from __future__ import annotations

import csv
from pathlib import Path


def _auc_wilcoxon(pos_scores: list[float], neg_scores: list[float]) -> float:
    """Compute AUROC via the Wilcoxon-Mann-Whitney statistic (O(n*m) but n is small)."""
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return 0.5
    concordant = sum(
        1 for p in pos_scores for n in neg_scores if p > n
    ) + 0.5 * sum(
        1 for p in pos_scores for n in neg_scores if p == n
    )
    return concordant / (n_pos * n_neg)


def _recall_at_k(labels: list[int], k: int) -> float:
    """Fraction of true positives in the top-k ranked items."""
    n_pos = sum(labels)
    if n_pos == 0:
        return 0.0
    top_k_pos = sum(labels[:k])
    return top_k_pos / n_pos


def run_retrospective_benchmark(
    amp_csv: str | Path,
    decoy_csv: str | Path,
    config_path: str | Path = "configs/pipeline.yaml",
    recall_ks: list[int] | None = None,
) -> dict:
    """Score known AMPs and shuffled decoys and compute AUROC + recall@k.

    Returns a dict with AUROC, per-k recall, and the full ranked list.
    """
    from openamp_foundry.features.physchem import compute_features
    from openamp_foundry.scoring.activity import activity_likeness_score
    from openamp_foundry.scoring.boman import boman_activity_score, gravy_score
    from openamp_foundry.scoring.ensemble import ensemble_score
    from openamp_foundry.scoring.novelty import novelty_score
    from openamp_foundry.scoring.safety import safety_score
    from openamp_foundry.scoring.synthesis import synthesis_feasibility_score
    from openamp_foundry.config import load_config

    config = load_config(config_path)
    weights = config["weights"]

    rows = []

    for path, true_label in [(amp_csv, 1), (decoy_csv, 0)]:
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                seq = row["sequence"].strip().upper()
                seq_id = row["id"]
                features = compute_features(seq)
                act = activity_likeness_score(features)
                safe = safety_score(features)
                synth = synthesis_feasibility_score(features, valid_sequence=True)
                nov, _ = novelty_score(seq, [])
                boman_act = boman_activity_score(seq)
                raw_scores = {
                    "activity": act, "safety": safe,
                    "synthesis": synth, "novelty": nov,
                    "boman_activity": boman_act,
                    "disagreement": abs(act - boman_act),
                }
                raw_scores["ensemble"] = ensemble_score(raw_scores, weights)
                rows.append({
                    "id": seq_id,
                    "sequence": seq,
                    "label": true_label,
                    "ensemble": raw_scores["ensemble"],
                    "activity": act,
                    "safety": safe,
                    "boman_activity": boman_act,
                    "hydrophobic_moment": features.get("hydrophobic_moment", 0.0),
                })

    rows.sort(key=lambda r: r["ensemble"], reverse=True)

    pos_scores = [r["ensemble"] for r in rows if r["label"] == 1]
    neg_scores = [r["ensemble"] for r in rows if r["label"] == 0]
    auroc = round(_auc_wilcoxon(pos_scores, neg_scores), 4)

    n_total = len(rows)
    n_pos = sum(r["label"] for r in rows)
    labels_ranked = [r["label"] for r in rows]

    if recall_ks is None:
        recall_ks = [10, 20, 44]
    recall = {f"recall_at_{k}": round(_recall_at_k(labels_ranked, k), 4) for k in recall_ks}

    random_auroc = 0.5
    interpretation = (
        "STRONG — model has meaningful discriminative power (AUROC > 0.70)"
        if auroc >= 0.70
        else "WEAK — model has modest signal; proceed with caution (AUROC 0.55–0.70)"
        if auroc >= 0.55
        else "POOR — model is near-random (AUROC < 0.55); do NOT proceed to synthesis"
    )

    return {
        "benchmark": "retrospective_auroc",
        "n_positives": n_pos,
        "n_negatives": n_total - n_pos,
        "n_total": n_total,
        "auroc": auroc,
        "random_auroc": random_auroc,
        "auroc_above_random": round(auroc - random_auroc, 4),
        **recall,
        "interpretation": interpretation,
        "design_note": (
            "Negatives are amino-acid-composition-matched shuffled decoys (RNG seed=42). "
            "AUROC reflects discrimination by ORDER-DEPENDENT features only "
            "(primarily hydrophobic moment). Composition-based features "
            "(charge, hydrophobic fraction, Boman index, GRAVY) are identical "
            "for each AMP/decoy pair and do NOT contribute to discrimination."
        ),
        "known_blind_spots": [
            "Melittin-like bent-helix peptides: hemolytic character not captured "
            "by simple 1D hydrophobic moment (Habermann 1972).",
            "Proline-rich AMPs (PR-39): activity relies on intracellular targets, "
            "not membrane disruption; hydrophobic moment is low but activity is real.",
        ],
        "top_ranked": rows[:10],
        "disclaimer": (
            "AUROC > 0.70 does NOT imply the nominated candidates are antimicrobial. "
            "It implies the model has some discriminative power over composition-matched "
            "controls. Wet-lab validation remains mandatory."
        ),
    }
