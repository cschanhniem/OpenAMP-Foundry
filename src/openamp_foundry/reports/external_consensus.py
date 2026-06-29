"""External predictor consensus aggregation.

Reads a CSV of per-candidate binary results from 5 external tools
(CAMPR4, AMPScanner v2, dbAMP 2.0, AntiCP 2.0, Macrel) and produces
a consensus verdict per candidate.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TOOLS = ["CAMPR4", "AMPScanner", "dbAMP", "AntiCP2", "Macrel"]
TOOLS_SHORT = ["CAMP", "AMPscn", "dbAMP", "AntiCP", "Macrel"]

# Safety-only predictors (hemolysis, toxicity) — tracked but not counted as AMP votes.
SAFETY_TOOLS = ["HAPPENN"]
SAFETY_TOOLS_SHORT = ["HAPPN"]


@dataclass
class ConsensusResult:
    candidate_id: str
    sequence: str
    votes: dict[str, bool]  # tool_name -> True (AMP) / False (non-AMP)
    n_positive: int
    n_tools: int
    consensus: str  # CONFIDENT / UNCERTAIN / WEAK
    mechanism_note: str


def compute_consensus(
    pilot_csv: str | Path, results_csv: str | Path
) -> list[ConsensusResult]:
    """Compute per-candidate consensus from filled-in predictor results.

    Args:
        pilot_csv: Path to pilot_panel.csv (must have 'candidate_id', 'sequence')
        results_csv: Path to external_predict_results.csv (Y/N per tool per candidate)

    Returns:
        List of ConsensusResult objects.
    """
    # Load pilot panel for sequences
    panel: dict[str, str] = {}
    with open(pilot_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cid = row.get("candidate_id", "")
            seq = row.get("sequence", "")
            if cid:
                panel[cid] = seq

    # Load results
    results: list[ConsensusResult] = []
    with open(results_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cid = row.get("candidate_id", "")
            votes: dict[str, bool] = {}
            n_positive = 0
            n_tools = 0
            for tool in TOOLS:
                val = row.get(tool, "").strip().upper()
                if val == "Y":
                    votes[tool] = True
                    n_positive += 1
                    n_tools += 1
                elif val == "N":
                    votes[tool] = False
                    n_tools += 1
                # blank/skip = not submitted

            total = max(n_tools, 1)
            frac = n_positive / total
            if frac >= 0.6:
                cons = "CONFIDENT"
            elif frac >= 0.4:
                cons = "UNCERTAIN"
            else:
                cons = "WEAK"

            seq = panel.get(cid, "")
            # Mechanism note for AntiCP2 (ACP-not-AMP caveat)
            note = ""
            amp_tools = [t for t, v in votes.items() if v]
            if "AntiCP2" in amp_tools:
                note = (
                    "AntiCP2 predicts anticancer (ACP) activity, not AMP directly. "
                    "Count ACP-positive as indirect supporting evidence only."
                )

            results.append(ConsensusResult(
                candidate_id=cid,
                sequence=seq,
                votes=votes,
                n_positive=n_positive,
                n_tools=n_tools,
                consensus=cons,
                mechanism_note=note,
            ))

    return results


def write_consensus_report(results: list[ConsensusResult], out_path: str | Path) -> None:
    """Write Markdown consensus report."""
    from datetime import datetime, timezone

    n_confident = sum(1 for r in results if r.consensus == "CONFIDENT")
    n_uncertain = sum(1 for r in results if r.consensus == "UNCERTAIN")
    n_weak = sum(1 for r in results if r.consensus == "WEAK")

    # Per-tool stats
    tool_positive: dict[str, int] = {}
    tool_total: dict[str, int] = {}
    for r in results:
        for tool, val in r.votes.items():
            tool_total[tool] = tool_total.get(tool, 0) + 1
            if val:
                tool_positive[tool] = tool_positive.get(tool, 0) + 1

    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    lines = [
        "# External Predictor Consensus Report",
        "",
        f"> **Generated:** {now}",
        f"> **Panel:** {len(results)} candidates",
        f"> **Tools:** {', '.join(TOOLS)}",
        "> **Consensus rule:** ≥60% positive = CONFIDENT, 40-59% = UNCERTAIN, <40% = WEAK",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Verdict | Count | Action |",
        "|---------|:-----:|--------|",
        f"| **CONFIDENT** (≥3/5) | {n_confident} | Proceed to synthesis |",
        f"| **UNCERTAIN** (2/5) | {n_uncertain} | Expert review required |",
        f"| **WEAK** (≤1/5) | {n_weak} | Do not synthesise without strong internal score |",
        "",
        "## Per-Tool Positive Rate",
        "",
        "| Tool | Positive | Total | Rate |",
        "|------|:--------:|:-----:|:----:|",
    ]
    for tool, short in zip(TOOLS + SAFETY_TOOLS, TOOLS_SHORT + SAFETY_TOOLS_SHORT):
        pos = tool_positive.get(tool, 0)
        tot = tool_total.get(tool, 0)
        tag = "(AMP)" if tool in TOOLS else "(safety)"
        rate = f"{pos}/{tot} ({100*pos/tot:.0f}%)" if tot else "N/A"
        lines.append(f"| {short} {tag} | {pos} | {tot} | {rate} |")

    lines.extend([
        "",
        "## Per-Candidate Results",
        "",
        "| Candidate | Sequence | CAMP | AMPscn | dbAMP | AntiCP | Macrel | HAPPN‡ | Agree | Verdict |",
        "|-----------|----------|:----:|:------:|:-----:|:------:|:------:|:-----:|:-----:|:--------:|",
    ])
    # Safety tool marks are shown in a separate column (not counted in consensus)
    all_col_tools = TOOLS + SAFETY_TOOLS
    for r in results:
        tool_marks = "".join(
            "Y" if r.votes.get(t, False) else ("N" if t in r.votes else ".")
            for t in all_col_tools
        )
        agree = f"{r.n_positive}/{r.n_tools}"
        seq_short = r.sequence[:25] + ("..." if len(r.sequence) > 25 else "")
        icon = {"CONFIDENT": "✅", "UNCERTAIN": "⚠️", "WEAK": "❌"}.get(r.consensus, "❓")
        lines.append(
            f"| {r.candidate_id} | {seq_short} | {tool_marks[0]} | {tool_marks[1]} | "
            f"{tool_marks[2]} | {tool_marks[3]} | {tool_marks[4]} | {tool_marks[5]} | "
            f"{agree} | {icon} {r.consensus} |"
        )

    lines.extend([
        "",
        "## Interpretation",
        "",
        "### CONFIDENT candidates",
        "At least 3 of 5 external tools predict antimicrobial activity. "
        "These have the highest synthesis priority.",
        "",
        "### UNCERTAIN candidates",
        "Only 2 of 5 tools agree. Expert review is required before synthesis. "
        "Consider the mechanistic basis for disagreement (e.g. AntiCP2 predicts "
        "anticancer peptides, not AMPs directly).",
        "",
        "### WEAK candidates",
        "0-1 of 5 tools predict antimicrobial activity. These should not be "
        "synthesised unless the internal pipeline score is exceptionally strong "
        "(ensemble > 0.85) with a clear mechanistic justification.",
        "",
        "## Caveats",
        "",
        "1. External tools are also computational — not wet-lab evidence.",
        "2. AntiCP 2.0 predicts **anticancer peptides (ACPs)**, not AMPs directly. "
        "ACP and AMP activity correlate but are not identical. Count ACP-positive "
        "as indirect supporting evidence only.",
        "3. Macrel has a known ONNX bug in local install (PR #77). Use the web server. "
        "(`big-data-biology.org/software/macrel`)",
        "4. **‡HAPPENN** predicts hemolysis risk, not AMP activity. HAPPENN column is shown "
        "but NOT counted in the consensus vote tally. HAPPENN-positive means increased "
        "hemolysis risk, which may conflict with the safety gate.",
        "5. Results not yet available for this panel — results table currently shows "
        "placeholder data. See `outputs/external_predict_checklist.md` for submission guide.",
    ])

    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def consensus_report_to_dict(results: list[ConsensusResult]) -> dict[str, Any]:
    """Convert consensus results to JSON-serialisable dict."""
    return {
        "status": "ok",
        "n_candidates": len(results),
        "n_confident": sum(1 for r in results if r.consensus == "CONFIDENT"),
        "n_uncertain": sum(1 for r in results if r.consensus == "UNCERTAIN"),
        "n_weak": sum(1 for r in results if r.consensus == "WEAK"),
        "results": [
            {
                "candidate_id": r.candidate_id,
                "consensus": r.consensus,
                "n_positive": r.n_positive,
                "n_tools": r.n_tools,
                "votes": {t: r.votes.get(t, None) for t in TOOLS},
            }
            for r in results
        ],
    }
