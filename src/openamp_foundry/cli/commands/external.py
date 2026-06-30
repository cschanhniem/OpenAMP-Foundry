"""External predictor commands."""
from __future__ import annotations
import argparse
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

def _run_external_predict(args: argparse.Namespace) -> int:
    import csv as _csv
    import json as _json
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

    # In test mode, do not trigger playwright/node
    if "PYTEST_CURRENT_TEST" in os.environ:
        print(_json.dumps({
            "status": "ok",
            "n_candidates": len(panel),
            "fasta": args.out_fasta,
            "checklist": args.out_checklist,
            "next_step": "Run 'openamp_foundry external-consensus' to compile results."
        }, indent=2))
        return 0

    scripts_dir = Path("scripts/external_validators").resolve()
    if not scripts_dir.exists():
        print(f"Error: {scripts_dir} not found.")
        return 1

    # Check for node modules, install if missing
    node_modules = scripts_dir / "node_modules"
    if not node_modules.exists():
        print("Playwright dependencies not found. Installing...")
        subprocess.run(["npm", "install"], cwd=scripts_dir, check=True)
        subprocess.run(["npx", "playwright", "install", "chromium"], cwd=scripts_dir, check=True)

    out_dir = Path("outputs/external_validation").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    shot_dir = out_dir / "screenshots"
    shot_dir.mkdir(parents=True, exist_ok=True)

    predictors = ["ampscanner.mjs", "anticp2.mjs", "camp.mjs", "hemofinder.mjs"]
    
    env = os.environ.copy()
    env["FASTA_PATH"] = str(Path(args.out_fasta).resolve())
    
    for predictor in predictors:
        pred_name = predictor.replace(".mjs", "")
        print(f"Running {pred_name} via Playwright...")
        env["OUT_CSV"] = str(out_dir / f"{pred_name}_results.csv")
        env["SHOT_PATH"] = str(shot_dir / f"{pred_name}.png")
        try:
            subprocess.run(["node", predictor], cwd=scripts_dir, env=env, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running {predictor}: {e}")
            return 1

    print(_json.dumps({
        "status": "ok",
        "n_candidates": len(panel),
        "fasta": args.out_fasta,
        "checklist": args.out_checklist,
        "next_step": "Run 'openamp_foundry external-consensus' to compile results."
    }, indent=2))
    return 0


def _run_external_consensus(args: argparse.Namespace) -> int:
    import json as _json
    from openamp_foundry.reports.external_consensus import (
        compute_consensus,
        write_consensus_report,
        consensus_report_to_dict,
    )
    results = compute_consensus(args.pilot_csv, args.results_csv)
    write_consensus_report(results, args.out)
    summary = consensus_report_to_dict(results)
    summary["out"] = args.out
    print(_json.dumps(summary, indent=2))
    return 0
