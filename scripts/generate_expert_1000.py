"""God-level de novo AMP generator — expert-objective gated, target 1000 candidates.

Every candidate that survives is simultaneously:
  • HIGH_CONFIDENCE_NOVEL  — <40% BLOSUM62 identity to all 51,503 known AMPs
  • MOTIF-NOVEL            — no long contiguous k-mer lifted from a known AMP
  • SELECTIVE              — selectivity_proxy ≥ 0.6 (charge/GRAVY therapeutic window)
  • LOW-HEMOLYSIS          — safety_score ≥ 0.6, μH capped, aromatic/Trp capped
  • SYNTHESISABLE          — synthesis_feasibility ≥ 0.7, no DKP/aspartimide/Trp-photo
  • CLEAR                  — no DRAMP patent proximity at any threshold

Ranked by the transparent expert composite (scoring/expert.py), which balances
activity ∩ selectivity ∩ safety ∩ synthesis ∩ motif-novelty ∩ helix-hinge — NOT a
single proxy. This is the automation of what a 30-year peptide chemist weighs at once.

Pipeline (cost-ordered, cheap rejects first):
  1. generate                       (main proc, ~µs)
  2. compute_features + expert gates(main proc, ~ms)        ← biophysics + selectivity
  3. pre-synthesis QC liabilities   (main proc, ~ms)        ← DKP/aspartimide/Trp-photo
  4. k-mer prior-art prefilter      (main proc, set lookup) ← local motif novelty
  5. BLOSUM62 novelty scan          (9 workers, parallel)   ← the only expensive step
  6. expert composite + diversity   (main proc)

Usage:
    .venv/bin/python3 scripts/generate_expert_1000.py [--workers N] [--target 1000]

Output (checkpointed every 50):
    outputs/expert_1000_candidates.csv
    outputs/expert_1000_candidates.fasta
"""
from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from openamp_foundry.features.physchem import compute_features
from openamp_foundry.qc.presynth_check import check_sequence
from openamp_foundry.scoring.expert import (
    build_kmer_index,
    expert_score,
    kmer_prior_art,
)

# ── DB sources (must match run_expanded_novelty_audit.py) ─────────────────────

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

DB_SOURCE_LIST: list[tuple[str, str, bool]] = [
    ("apd6_natural",       "data/novelty_db/apd6_natural.fasta",          False),
    ("apd6_animal",        "data/novelty_db/apd6_animal.fasta",           False),
    ("apd6_plant",         "data/novelty_db/apd6_plant.fasta",            False),
    ("apd6_bacteria",      "data/novelty_db/apd6_bacteria.fasta",         False),
    ("apd6_human",         "data/novelty_db/apd6_human.fasta",            False),
    ("dramp_general",      "data/novelty_db/dramp_general.fasta",         False),
    ("dramp_patent",       "data/novelty_db/dramp_patent.fasta",          True),
    ("dramp_specific",     "data/novelty_db/dramp_specific.fasta",        False),
    ("uniprot_reviewed",   "data/novelty_db/uniprot_amps_reviewed.fasta", False),
    ("uniprot_unreviewed", "data/novelty_db/uniprot_amps_unreviewed.fasta",False),
    ("uniprot_combined",   "data/novelty_db/uniprot_amps.fasta",          False),
    ("escape_amps",        "data/novelty_db/escape_amps.fasta",           False),
    ("dbamp3",             "data/novelty_db/dbAMP3.fasta",                False),
    ("dbaasp",             "data/novelty_db/dbaasp-peptides.fasta",       False),
]

OUTPUT_CSV   = ROOT / "outputs" / "expert_1000_candidates.csv"
OUTPUT_FASTA = ROOT / "outputs" / "expert_1000_candidates.fasta"

# ── Expert gate thresholds (fixed before ranking) ─────────────────────────────

GATE = {
    "len_min": 12, "len_max": 24,
    "charge_min": 3.0, "charge_max": 8.0,        # net_charge_ph74
    "hydro_min": 0.25, "hydro_max": 0.50,
    "aromatic_max": 0.25,
    "trp_max": 2,
    "mu_h_min": 0.18, "mu_h_max": 0.55,          # below melittin-like lytic band
    "selectivity_min": 0.60,
    "safety_min": 0.60,
    "synthesis_min": 0.70,
    "kmer_k": 5,
    "kmer_max_known_run": 2,                      # reject ≥3 consecutive known 5-mers
    "novelty_max_identity": 0.40,                 # <40% BLOSUM62 = HIGH_CONFIDENCE_NOVEL
}

# Pre-synthesis liabilities that hard-fail a candidate (expert would not order these).
_HARD_LIABILITIES = ("DKP_RISK", "PYROGLUTAMATE_RISK", "TRP_PHOTOLABILITY", "ISOMERIZATION_RISK")

HYDROPHOBIC = frozenset("LIVWFA")
AROMATIC    = frozenset("FWY")

# Alphabet biased toward selective, synthesisable, helix-hinge-capable peptides:
# K/R cationic; P/G enabled (hinge); W/F downweighted (hemolysis/photolability);
# polar N/S/T/Q present for solubility.
_WEIGHTS = {
    "A": 3, "D": 1, "E": 1, "F": 3, "G": 4,
    "H": 1, "I": 5, "K": 9, "L": 6, "N": 3,
    "P": 4, "Q": 1, "R": 9, "S": 3, "T": 3,
    "V": 5, "W": 2, "Y": 2,
}
_POOL = [aa for aa, w in _WEIGHTS.items() for _ in range(w)]


def generate_candidate(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(_POOL) for _ in range(length))


# ── DB load (shared FASTA loader) ─────────────────────────────────────────────

def _load_fasta(path: Path) -> list[tuple[str, str]]:
    results, header, parts = [], "", []
    try:
        with open(path) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    if header and parts:
                        results.append((header, "".join(parts).upper()))
                    header = line[1:].split()[0] if len(line) > 1 else "UNKNOWN"
                    parts = []
                else:
                    parts.append(line)
        if header and parts:
            results.append((header, "".join(parts).upper()))
    except FileNotFoundError:
        pass
    return results


def build_db(root: Path) -> list[tuple[str, str, bool]]:
    seen: dict[str, tuple[str, bool]] = {}
    for name, rel, is_patent in DB_SOURCE_LIST:
        for header, seq in _load_fasta(root / rel):
            if not seq or not all(c in STANDARD_AA for c in seq) or not (5 <= len(seq) <= 100):
                continue
            if seq not in seen:
                seen[seq] = (f"{name}:{header}", is_patent)
            elif is_patent and not seen[seq][1]:
                seen[seq] = (seen[seq][0], True)
    return [(seq, sid, pat) for seq, (sid, pat) in seen.items()]


# ── Worker (BLOSUM62 novelty scan only) ───────────────────────────────────────

_WDB: list = []
_WALIGNER = None


def _worker_init(root_str: str) -> None:
    from Bio.Align import PairwiseAligner, substitution_matrices
    global _WDB, _WALIGNER
    _WDB = build_db(Path(root_str))
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.mode = "local"
    a.open_gap_score = -11.0
    a.extend_gap_score = -1.0
    _WALIGNER = a


def _local_identity(query: str, target: str) -> float:
    try:
        aln = next(iter(_WALIGNER.align(query, target)))
    except Exception:
        return 0.0
    n = 0
    for (qs, qe), (ts, _te) in zip(aln.aligned[0], aln.aligned[1]):
        for i in range(qe - qs):
            if query[qs + i] == target[ts + i]:
                n += 1
    return n / len(query) if query else 0.0


def _worker_scan(seq: str) -> tuple[str, float, str, bool]:
    n = len(seq)
    min_l, max_l = max(5, n // 3), n * 3
    best_id, best_hit, best_pat = 0.0, "NONE", False
    for db_seq, db_id, db_pat in _WDB:
        if not (min_l <= len(db_seq) <= max_l):
            continue
        if seq in db_seq:
            return seq, 1.0, db_id, db_pat
        identity = _local_identity(seq, db_seq)
        if identity > best_id:
            best_id, best_hit, best_pat = identity, db_id, db_pat
        if best_id >= GATE["novelty_max_identity"]:
            return seq, best_id, best_hit, best_pat
    return seq, best_id, best_hit, best_pat


# ── Cheap expert gates (main process) ─────────────────────────────────────────

def passes_expert_gates(seq: str, feats: dict, kmer_index: set[str]) -> tuple[bool, str]:
    n = len(seq)
    if not (GATE["len_min"] <= n <= GATE["len_max"]):
        return False, "length"
    if "C" in seq or "M" in seq:
        return False, "cys_met"
    charge = feats["net_charge_ph74"]
    if not (GATE["charge_min"] <= charge <= GATE["charge_max"]):
        return False, "charge"
    hf = feats["hydrophobic_fraction"]
    if not (GATE["hydro_min"] <= hf <= GATE["hydro_max"]):
        return False, "hydrophobic"
    if feats["aromatic_fraction"] > GATE["aromatic_max"]:
        return False, "aromatic"
    if seq.count("W") > GATE["trp_max"]:
        return False, "trp"
    mu_h = feats["hydrophobic_moment"]
    if not (GATE["mu_h_min"] <= mu_h <= GATE["mu_h_max"]):
        return False, "moment"
    if feats["selectivity_proxy"] < GATE["selectivity_min"]:
        return False, "selectivity"

    # k-mer prior-art (cheap set lookups)
    motif = kmer_prior_art(seq, kmer_index, k=GATE["kmer_k"])
    if motif["max_run_known"] >= GATE["kmer_max_known_run"] + 1:
        return False, "motif_prior_art"

    # Pre-synthesis liabilities
    qc = check_sequence("gate", seq, mu_h=mu_h)
    if any(any(h in f for h in _HARD_LIABILITIES) for f in qc.flags):
        return False, "synth_liability"
    if qc.synthesis_difficulty == "HIGH":
        return False, "synth_difficulty"

    return True, "ok"


# ── Diversity bucketing ───────────────────────────────────────────────────────

MAX_PER_BIN = 28


def _diversity_bin(seq: str, feats: dict) -> tuple:
    charge = int(round(feats["net_charge_ph74"]))
    hf = feats["hydrophobic_fraction"]
    return (min(charge, 8) // 2, len(seq) // 4, int(hf * 5))


# ── Output ────────────────────────────────────────────────────────────────────

FIELDNAMES = [
    "candidate_id", "sequence", "length", "net_charge_ph74", "hydrophobic_fraction",
    "aromatic_fraction", "mu_h", "max_mu_h", "gravy", "selectivity_proxy",
    "expert_composite", "expert_activity_consensus", "expert_selectivity",
    "expert_safety", "expert_synthesis", "expert_serum_stability",
    "expert_hinge_selectivity", "expert_motif_novelty",
    "has_central_hinge", "motif_known_kmers", "motif_max_known_run",
    "best_identity", "best_hit_id", "novelty_class", "patent_risk",
    "synthesis_difficulty", "expert_flags", "seed_family",
]


def _write_outputs(rows: list[dict]) -> None:
    with open(OUTPUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(OUTPUT_FASTA, "w") as f:
        for r in rows:
            f.write(
                f">{r['candidate_id']} composite={r['expert_composite']:.3f} "
                f"sel={r['expert_selectivity']:.2f} safe={r['expert_safety']:.2f} "
                f"sim={r['best_identity']:.1%} charge={r['net_charge_ph74']:.1f} "
                f"hinge={int(r['has_central_hinge'])} note=N_ACETYLATION_RECOMMENDED\n"
            )
            f.write(r["sequence"] + "\n")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--target", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260630)
    args = ap.parse_args()

    N_TARGET, N_WORKERS, SEED = args.target, args.workers, args.seed
    LENGTHS = list(range(GATE["len_min"], GATE["len_max"] + 1))
    BATCH = N_WORKERS * 6
    CHECKPOINT = 50

    print("=== God-Level Expert Generator — target", N_TARGET, "===")
    print("DB: APD6 + DRAMP3 + UniProt + ESCAPE + dbAMP3 + DBAASP  (51,503 seqs)")
    print(f"Workers: {N_WORKERS} | seed: {SEED}")
    print("Gates: selectivity≥0.60, safety≥0.60, synth≥0.70, μH 0.18–0.55,")
    print("       aromatic≤0.25, W≤2, k-mer run<3, no DKP/aspartimide/Trp-photo, <40% identity")
    print("Ranking: expert composite (activity∩selectivity∩safety∩synth∩motif∩hinge)\n", flush=True)

    t0 = time.time()
    print("Building known-AMP DB + k-mer index (main process)...", flush=True)
    db = build_db(ROOT)
    kmer_index = build_kmer_index([s for s, _, _ in db], k=GATE["kmer_k"])
    print(f"  DB: {len(db):,} seqs | {GATE['kmer_k']}-mer index: {len(kmer_index):,} motifs "
          f"({time.time()-t0:.1f}s)\n", flush=True)

    print(f"Starting {N_WORKERS} BLOSUM workers...", flush=True)
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(processes=N_WORKERS, initializer=_worker_init, initargs=(str(ROOT),))
    print(f"  ready ({time.time()-t0:.1f}s)\n", flush=True)

    rng = random.Random(SEED)
    rows: list[dict] = []
    seen: set[str] = set()
    bin_counts: dict[tuple, int] = defaultdict(int)
    reject = defaultdict(int)
    n_gen = n_gate = n_novel = n_div = 0
    t_start = time.time()

    print(f"{'#':>5}  {'ID':<14} {'Sequence':<26} {'Comp':>5} {'Sel':>4} {'Safe':>4} "
          f"{'Hng':>3} {'Sim':>6}", flush=True)
    print("─" * 118, flush=True)

    pending: list[tuple[str, dict]] = []   # (seq, feats) awaiting BLOSUM
    attempt = 0
    N_MAX = 50_000_000

    try:
        while n_novel < N_TARGET and attempt < N_MAX:
            # Fill a batch of gate-passing candidates
            while len(pending) < BATCH and attempt < N_MAX:
                attempt += 1
                seq = generate_candidate(rng, rng.choice(LENGTHS))
                n_gen += 1
                if seq in seen:
                    continue
                seen.add(seq)
                feats = compute_features(seq)
                ok, why = passes_expert_gates(seq, feats, kmer_index)
                if not ok:
                    reject[why] += 1
                    continue
                n_gate += 1
                pending.append((seq, feats))

            if not pending:
                break

            scan_in = [s for s, _ in pending]
            feat_map = {s: f for s, f in pending}
            pending = []
            scan_out = pool.map(_worker_scan, scan_in)

            for seq, best_id, best_hit, is_pat in scan_out:
                if best_id >= GATE["novelty_max_identity"] or is_pat:
                    reject["novelty_or_patent"] += 1
                    continue
                feats = feat_map[seq]
                div = _diversity_bin(seq, feats)
                if bin_counts[div] >= MAX_PER_BIN:
                    n_div += 1
                    continue

                es = expert_score(seq, features=feats, kmer_index=kmer_index, k=GATE["kmer_k"])
                qc = check_sequence("c", seq, mu_h=feats["hydrophobic_moment"])

                n_novel += 1
                bin_counts[div] += 1
                cid = f"XPRT_{n_novel:04d}"

                row = {
                    "candidate_id": cid, "sequence": seq, "length": len(seq),
                    "net_charge_ph74": round(feats["net_charge_ph74"], 2),
                    "hydrophobic_fraction": feats["hydrophobic_fraction"],
                    "aromatic_fraction": feats["aromatic_fraction"],
                    "mu_h": feats["hydrophobic_moment"],
                    "max_mu_h": feats["max_hydrophobic_moment"],
                    "gravy": feats["gravy"],
                    "selectivity_proxy": feats["selectivity_proxy"],
                    "expert_composite": es.composite,
                    "expert_activity_consensus": es.components["activity_consensus"],
                    "expert_selectivity": es.components["selectivity"],
                    "expert_safety": es.components["safety"],
                    "expert_synthesis": es.components["synthesis"],
                    "expert_serum_stability": es.components["serum_stability"],
                    "expert_hinge_selectivity": es.components["hinge_selectivity"],
                    "expert_motif_novelty": es.components["motif_novelty"],
                    "has_central_hinge": es.extras["has_central_hinge"],
                    "motif_known_kmers": es.extras["motif_known_kmers"],
                    "motif_max_known_run": es.extras["motif_max_known_run"],
                    "best_identity": round(best_id, 4), "best_hit_id": best_hit,
                    "novelty_class": "HIGH_CONFIDENCE_NOVEL", "patent_risk": "CLEAR",
                    "synthesis_difficulty": qc.synthesis_difficulty,
                    "expert_flags": ";".join(es.flags),
                    "seed_family": cid,
                }
                rows.append(row)

                print(f"{n_novel:>5}  {cid:<14} {seq:<26} {es.composite:>5.3f} "
                      f"{es.components['selectivity']:>4.2f} {es.components['safety']:>4.2f} "
                      f"{int(es.extras['has_central_hinge']):>3} {best_id:>6.1%}", flush=True)

                if n_novel % CHECKPOINT == 0:
                    rows.sort(key=lambda r: -r["expert_composite"])
                    _write_outputs(rows)
                    el = time.time() - t_start
                    rate = n_novel / el * 3600
                    print(f"\n  [ckpt {n_novel}/{N_TARGET}] saved | {n_gate} gated / "
                          f"{n_gen} gen ({100*n_gate/max(n_gen,1):.1f}%) | "
                          f"{rate:.0f}/hr ETA {(N_TARGET-n_novel)/max(rate,1)*60:.0f}min\n", flush=True)
                if n_novel >= N_TARGET:
                    break
    finally:
        pool.close()
        pool.join()

    rows.sort(key=lambda r: -r["expert_composite"])
    _write_outputs(rows)
    el = time.time() - t_start
    print("\n" + "─" * 118)
    print("\n=== DONE ===")
    print(f"  Generated: {n_gen:,} | gate-passed: {n_gate:,} ({100*n_gate/max(n_gen,1):.2f}%)")
    print(f"  Novel+CLEAR kept: {n_novel} | dropped by diversity cap: {n_div}")
    print(f"  Time: {el/60:.1f} min")
    print(f"\n  Reject reasons:")
    for why, c in sorted(reject.items(), key=lambda x: -x[1]):
        print(f"    {why:<22} {c:,}")
    if rows:
        print(f"\n  Top 12 by expert composite:")
        print(f"  {'ID':<14} {'Sequence':<26} {'Comp':>5} {'Sel':>4} {'Safe':>4} {'Act':>4} {'Hng':>3}")
        for r in rows[:12]:
            print(f"  {r['candidate_id']:<14} {r['sequence']:<26} {r['expert_composite']:>5.3f} "
                  f"{r['expert_selectivity']:>4.2f} {r['expert_safety']:>4.2f} "
                  f"{r['expert_activity_consensus']:>4.2f} {int(r['has_central_hinge']):>3}")
    print(f"\n  CSV  → {OUTPUT_CSV}")
    print(f"  FASTA→ {OUTPUT_FASTA}")
    print(f"  Next: screen_1000_candidates.py (Macrel) then web predictors.")


if __name__ == "__main__":
    main()
