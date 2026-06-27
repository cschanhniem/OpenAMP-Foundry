from __future__ import annotations

import argparse
import json
from pathlib import Path

from openamp_foundry.evidence.schemas import validate_json_schema
from openamp_foundry.pipeline import run_ranking_pipeline
from openamp_foundry.utils.io import read_json


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="openamp-foundry")
    sub = parser.add_subparsers(dest="command", required=True)

    rank = sub.add_parser("rank", help="Rank candidate peptides and generate evidence.")
    rank.add_argument("--candidates", required=True)
    rank.add_argument("--references", required=False)
    rank.add_argument("--out", required=True)
    rank.add_argument("--report", required=False)
    rank.add_argument("--cert-dir", required=False)
    rank.add_argument("--manifest", required=False)
    rank.add_argument("--config", default="configs/pipeline.yaml")

    validate = sub.add_parser("validate", help="Validate a candidate certificate against JSON schema.")
    validate.add_argument("--certificate", required=True)
    validate.add_argument("--schema", required=True)

    bench = sub.add_parser("bench", help="Run benchmark and leakage checks.")
    bench_sub = bench.add_subparsers(dest="bench_command", required=True)

    leakage = bench_sub.add_parser("leakage", help="Find near-duplicate candidates in references.")
    leakage.add_argument("--candidates", required=True)
    leakage.add_argument("--references", required=True)
    leakage.add_argument("--threshold", type=float, default=0.90)
    leakage.add_argument("--out", required=False, help="Optional JSON output path.")

    baseline = bench_sub.add_parser(
        "baseline",
        help="Evaluate pipeline recall vs random baseline on a labelled set.",
    )
    baseline.add_argument("--candidates", required=True, help="CSV of all candidates to score.")
    baseline.add_argument("--references", required=False, help="Reference CSV for novelty scoring.")
    baseline.add_argument(
        "--positives",
        required=True,
        help="CSV of known-positive (active) peptides. IDs must match candidates CSV.",
    )
    baseline.add_argument(
        "--k",
        type=int,
        nargs="+",
        default=None,
        help="Recall@k cutoffs (default: auto from dataset size).",
    )
    baseline.add_argument("--config", default="configs/pipeline.yaml")
    baseline.add_argument("--out", required=False, help="Optional JSON output path.")

    generate = sub.add_parser(
        "generate-batch",
        help=(
            "Generate a candidate pool by conservative mutation of seed sequences. "
            "Output is a CSV suitable for the 'rank' command. "
            "This is a toy exploration tool — no biological activity is implied."
        ),
    )
    generate.add_argument(
        "--seeds",
        required=True,
        help="CSV of seed sequences (columns: id, sequence, source)",
    )
    generate.add_argument(
        "--out",
        required=True,
        help="Output CSV path for the candidate pool.",
    )
    generate.add_argument(
        "--n-double",
        type=int,
        default=25,
        help="Double-substitution variants per seed (default: 25)",
    )
    generate.add_argument(
        "--n-charge",
        type=int,
        default=12,
        help="Charge-enhanced variants per seed (default: 12)",
    )
    generate.add_argument(
        "--rng-seed",
        type=int,
        default=2024,
        help="RNG seed for reproducibility (default: 2024)",
    )

    pilot = sub.add_parser(
        "pilot-panel",
        help=(
            "Select a first-synthesis-wave pilot panel from a ranked JSONL file. "
            "Picks n candidates (default 20) maximising ensemble score, minimising scorer "
            "disagreement, and ensuring at least one representative per seed template."
        ),
    )
    pilot.add_argument(
        "--ranked",
        required=True,
        help="Ranked JSONL file (output of the 'rank' command).",
    )
    pilot.add_argument(
        "--n",
        type=int,
        default=20,
        help="Panel size (default: 20).",
    )
    pilot.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV path (synthesis-ready format).",
    )
    pilot.add_argument(
        "--out-md",
        required=False,
        help="Optional output path for human-readable markdown panel.",
    )

    validate_scoring = sub.add_parser(
        "validate-scoring",
        help=(
            "Retrospective AUROC benchmark: known AMPs vs background random peptides. "
            "AUROC > 0.70 = model passes Gate 1 (proceed to synthesis). "
            "Run before committing to wet-lab spend."
        ),
    )
    validate_scoring.add_argument(
        "--amp-csv",
        default="examples/validation/known_amps.csv",
        help="CSV of known AMPs with 'id' and 'sequence' columns (label=1).",
    )
    validate_scoring.add_argument(
        "--decoy-csv",
        default="examples/validation/random_background.csv",
        help=(
            "CSV of decoy peptides (label=0). "
            "Default: background-frequency random peptides (standard benchmark). "
            "Use examples/validation/scrambled_decoys.csv for the stricter "
            "composition-matched shuffle test."
        ),
    )
    validate_scoring.add_argument(
        "--benchmark-type",
        choices=["standard", "strict"],
        default="standard",
        help=(
            "standard: AMPs vs background random peptides (primary synthesis gate). "
            "strict: AMPs vs composition-matched shuffled decoys (order-sensitivity test)."
        ),
    )
    validate_scoring.add_argument("--config", default="configs/pipeline.yaml")
    validate_scoring.add_argument(
        "--out",
        required=False,
        help="Optional JSON output path.",
    )

    external_predict = sub.add_parser(
        "external-predict",
        help=(
            "Generate FASTA and submission checklist for external AMP prediction tools "
            "(CAMPR4, AMPScanner v2, dbAMP). Must be submitted manually — no API calls made."
        ),
    )
    external_predict.add_argument(
        "--pilot-csv",
        required=True,
        help="Pilot panel CSV (output of 'pilot-panel' command).",
    )
    external_predict.add_argument(
        "--out-fasta",
        default="outputs/pilot_panel.fasta",
        help="Output FASTA file for tool submission.",
    )
    external_predict.add_argument(
        "--out-checklist",
        default="outputs/external_predict_checklist.md",
        help="Output markdown checklist for recording results.",
    )

    pilot_confident = sub.add_parser(
        "pilot-confident",
        help=(
            "Filter a pilot panel to candidates confirmed by ≥2 external predictors. "
            "Provide the comma-separated IDs of confirmed candidates via --keep."
        ),
    )
    pilot_confident.add_argument(
        "--pilot-csv",
        required=True,
        help="Pilot panel CSV (output of 'pilot-panel' command).",
    )
    pilot_confident.add_argument(
        "--keep",
        required=True,
        help="Comma-separated candidate IDs to retain (from external predictor results).",
    )
    pilot_confident.add_argument(
        "--out",
        default="outputs/confident_panel",
        help="Output path prefix (will write .csv and .md).",
    )

    batch_pack = sub.add_parser(
        "batch-pack",
        help=(
            "Generate Phase 3 batch pack reports (diversity, novelty, toxicity, synthesis) "
            "from a ranked JSONL file produced by 'rank'."
        ),
    )
    batch_pack.add_argument(
        "--ranked",
        required=True,
        help="Ranked JSONL file (output of the 'rank' command).",
    )
    batch_pack.add_argument(
        "--out-json",
        required=True,
        help="Output path for machine-readable batch pack JSON.",
    )
    batch_pack.add_argument(
        "--out-md",
        required=False,
        help="Optional output path for human-readable markdown report.",
    )
    batch_pack.add_argument(
        "--diversity-threshold",
        type=float,
        default=0.80,
        help="Similarity threshold for diversity clustering (default: 0.80)",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "rank":
        run_ranking_pipeline(
            candidate_path=args.candidates,
            reference_path=args.references,
            out_path=args.out,
            report_path=args.report,
            cert_dir=args.cert_dir,
            config_path=args.config,
            manifest_path=args.manifest,
        )
        print(json.dumps({"status": "ok", "out": args.out, "report": args.report}, indent=2))
        return 0

    if args.command == "validate":
        payload = read_json(args.certificate)
        validate_json_schema(payload, Path(args.schema))
        print(json.dumps({"status": "valid", "certificate": args.certificate}, indent=2))
        return 0

    if args.command == "bench":
        return _run_bench(args)

    if args.command == "generate-batch":
        return _run_generate_batch(args)

    if args.command == "pilot-panel":
        return _run_pilot_panel(args)

    if args.command == "validate-scoring":
        return _run_validate_scoring(args)

    if args.command == "external-predict":
        return _run_external_predict(args)

    if args.command == "pilot-confident":
        return _run_pilot_confident(args)

    if args.command == "batch-pack":
        return _run_batch_pack(args)

    parser.error("unknown command")
    return 2


def _run_bench(args: argparse.Namespace) -> int:
    from openamp_foundry.data.loaders import load_candidates_csv
    from openamp_foundry.utils.io import write_json

    if args.bench_command == "leakage":
        from openamp_foundry.benchmark.leakage import find_near_duplicates

        candidates = load_candidates_csv(args.candidates)
        references = load_candidates_csv(args.references)
        hits = find_near_duplicates(candidates, references, threshold=args.threshold)
        result = {
            "status": "ok",
            "threshold": args.threshold,
            "near_duplicate_count": len(hits),
            "near_duplicates": hits,
            "warning": (
                "Near-duplicates detected. If these candidates were used for training or "
                "scoring baseline models, benchmark results may be inflated."
            ) if hits else None,
        }
        if args.out:
            write_json(args.out, result)
        print(json.dumps(result, indent=2))
        return 0

    if args.bench_command == "baseline":
        from openamp_foundry.benchmark.evaluate import benchmark_summary
        from openamp_foundry.pipeline import score_candidates

        scored, _ = score_candidates(
            candidate_path=args.candidates,
            reference_path=args.references,
            config_path=args.config,
        )
        positives = load_candidates_csv(args.positives)
        positive_ids = {p.candidate_id for p in positives}
        summary = benchmark_summary(scored, positive_ids, ks=args.k)
        result = {"status": "ok", **summary}
        if args.out:
            write_json(args.out, result)
        print(json.dumps(result, indent=2))
        return 0

    return 2


def _run_pilot_panel(args: argparse.Namespace) -> int:
    import json as _json
    from datetime import datetime, timezone

    from openamp_foundry.reports.pilot_panel import write_pilot_csv, write_pilot_markdown
    from openamp_foundry.selection.pilot import select_pilot_panel

    ranked_path = Path(args.ranked)
    if not ranked_path.exists():
        print(_json.dumps({"status": "error", "message": f"File not found: {args.ranked}"}))
        return 1

    candidates = []
    with open(ranked_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                row = _json.loads(line)
                if row.get("selected"):
                    candidates.append(row)

    panel = select_pilot_panel(candidates, n=args.n)
    generated_at = datetime.now(timezone.utc).isoformat()

    write_pilot_csv(panel, args.out_csv)
    if args.out_md:
        write_pilot_markdown(panel, args.out_md, generated_at=generated_at)

    seeds = sorted({c.get("seed", "") for c in panel})
    n_consensus = sum(
        1 for c in panel if c.get("scores", {}).get("disagreement", 1.0) < 0.20
    )
    print(_json.dumps({
        "status": "ok",
        "n_nominees": len(candidates),
        "n_panel": len(panel),
        "seeds_represented": seeds,
        "n_dual_scorer_consensus": n_consensus,
        "out_csv": args.out_csv,
        "out_md": args.out_md,
        "disclaimer": (
            "No antimicrobial activity has been demonstrated. "
            "Human expert review required before synthesis."
        ),
    }, indent=2))
    return 0


def _run_validate_scoring(args: argparse.Namespace) -> int:
    import json as _json
    from openamp_foundry.benchmark.retrospective import run_retrospective_benchmark
    from openamp_foundry.utils.io import write_json

    benchmark_type = getattr(args, "benchmark_type", "standard")
    result = run_retrospective_benchmark(
        amp_csv=args.amp_csv,
        decoy_csv=args.decoy_csv,
        config_path=args.config,
        benchmark_type=benchmark_type,
    )
    if args.out:
        write_json(args.out, result)
    summary = {
        "status": "ok",
        "auroc": result["auroc"],
        "auroc_above_random": result["auroc_above_random"],
        "recall_at_10": result.get("recall_at_10"),
        "recall_at_20": result.get("recall_at_20"),
        "recall_at_44": result.get("recall_at_44"),
        "interpretation": result["interpretation"],
        "out": args.out,
    }
    print(_json.dumps(summary, indent=2))
    return 0


def _run_external_predict(args: argparse.Namespace) -> int:
    import csv as _csv
    import json as _json
    from datetime import datetime, timezone
    from openamp_foundry.reports.external_predict import (
        write_external_predict_checklist,
        write_pilot_fasta,
    )

    panel = []
    with open(args.pilot_csv, newline="", encoding="utf-8") as f:
        for row in _csv.DictReader(f):
            panel.append(row)

    generated_at = datetime.now(timezone.utc).isoformat()
    write_pilot_fasta(panel, args.out_fasta)
    write_external_predict_checklist(
        panel, fasta_path=args.out_fasta, out_path=args.out_checklist,
        generated_at=generated_at,
    )
    print(_json.dumps({
        "status": "ok",
        "n_candidates": len(panel),
        "fasta": args.out_fasta,
        "checklist": args.out_checklist,
        "next_step": (
            f"Submit {args.out_fasta} to CAMPR4, AMPScanner v2, and dbAMP. "
            f"Fill in {args.out_checklist}. Then run 'make pilot-confident'."
        ),
    }, indent=2))
    return 0


def _run_pilot_confident(args: argparse.Namespace) -> int:
    import csv as _csv
    import json as _json
    from datetime import datetime, timezone
    from openamp_foundry.reports.external_predict import write_confident_panel

    panel = []
    with open(args.pilot_csv, newline="", encoding="utf-8") as f:
        for row in _csv.DictReader(f):
            panel.append(row)

    keep_ids = [cid.strip() for cid in args.keep.split(",") if cid.strip()]
    generated_at = datetime.now(timezone.utc).isoformat()
    confident = write_confident_panel(panel, keep_ids, out_path=args.out, generated_at=generated_at)

    print(_json.dumps({
        "status": "ok",
        "n_input": len(panel),
        "n_confident": len(confident),
        "out_csv": args.out + ".csv",
        "out_md": args.out + ".md",
        "disclaimer": "Confident candidates still require human expert review and biosafety sign-off.",
    }, indent=2))
    return 0


def _run_batch_pack(args: argparse.Namespace) -> int:
    from openamp_foundry.reports.batch_pack import generate_batch_pack, write_batch_pack_markdown
    from openamp_foundry.utils.io import write_json

    pack = generate_batch_pack(
        ranked_jsonl_path=args.ranked,
        diversity_threshold=args.diversity_threshold,
    )
    write_json(args.out_json, pack)
    if args.out_md:
        write_batch_pack_markdown(pack, args.out_md)

    print(json.dumps({
        "status": "ok",
        "n_selected": pack["summary"]["n_candidates_selected"],
        "n_clusters": pack["summary"]["n_diversity_clusters"],
        "mean_novelty": pack["summary"]["mean_novelty"],
        "mean_safety": pack["summary"]["mean_safety"],
        "mean_synthesis": pack["summary"]["mean_synthesis"],
        "out_json": args.out_json,
        "out_md": args.out_md,
    }, indent=2))
    return 0


def _run_generate_batch(args: argparse.Namespace) -> int:
    import csv

    from openamp_foundry.generators.template_mutator import generate_candidate_pool

    seeds_path = Path(args.seeds)
    if not seeds_path.exists():
        print(json.dumps({"status": "error", "message": f"Seeds file not found: {args.seeds}"}))
        return 1

    seed_ids: list[str] = []
    seed_seqs: list[str] = []
    with open(seeds_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            seed_ids.append(row["id"])
            seed_seqs.append(row["sequence"].strip().upper())

    pool = generate_candidate_pool(
        seed_sequences=seed_seqs,
        seed_ids=seed_ids,
        n_double=args.n_double,
        n_charge_enhance=args.n_charge,
        rng_seed=args.rng_seed,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["id", "sequence", "source"])
        writer.writeheader()
        writer.writerows(pool)

    print(json.dumps({
        "status": "ok",
        "n_seeds": len(seed_ids),
        "n_candidates_generated": len(pool),
        "out": str(out_path),
        "disclaimer": (
            "Generated candidates are toy conservative-substitution variants. "
            "They have no demonstrated biological activity. "
            "Run 'rank' to score and filter them."
        ),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
