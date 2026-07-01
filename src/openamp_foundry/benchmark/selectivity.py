"""Selectivity benchmark: can the pipeline distinguish selective AMPs from hemolytic AMPs?

This benchmark tests a fundamentally different question from the retrospective AUROC:
instead of "can we distinguish AMPs from random sequences?", it asks "can we distinguish
AMPs that are safe (low hemolysis) from AMPs that are dangerous (high hemolysis)?"

This is the critical benchmark for the safety-first objective. An AMP that kills bacteria
but also lyses RBCs is not a useful therapeutic candidate. If the pipeline cannot separate
selective from hemolytic AMPs, its safety filtering is not doing its job.

Reference panel: ``examples/validation/selectivity_panel.csv``
- 10 selective AMPs (literature-reported low hemolysis / HC50 > 100 uM or > 1 mg/mL)
- 8 hemolytic AMPs (literature-reported high hemolysis / HC50 < 50 uM)

Key metrics:
- AUROC of safety_score for selective-vs-hemolytic discrimination
- AUROC of selectivity_proxy for the same task
- AUROC of naive baselines (hydrophobic_fraction, hydrophobic_moment, charge_density)
- Per-peptide score table with flags for known blind spots

Literature basis for hemolysis classification:
- Melittin: HC50 ~ 5 ug/mL (Habermann 1972) — strongly hemolytic
- Mastoparan-X: hemolytic at >50 uM (Higashijima et al. 1990)
- BMAP-28: HC50 ~ 5 uM (Skerlavaj et al. 1996) — strongly hemolytic
- Magainin-2: HC50 > 1 mg/mL (~250 uM) (Zasloff 1987) — selective
- Cecropin-A: low hemolysis even at high concentration (Steiner et al. 1981)
- Buforin-II: non-lytic mechanism, no membrane disruption (Park et al. 2000)
- Histatin-5: low hemolysis, intracellular target (Oppenheim et al. 1988)
- Plectasin: low hemolysis, fungal defensin (Mygind et al. 2005)

All hemolysis classifications are literature-derived. The benchmark is computational only;
it tests whether physicochemical proxies correlate with literature hemolysis classifications.
It does not predict hemolysis for novel sequences.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class SelectivityRow:
    id: str
    sequence: str
    family: str
    hemolysis_class: str  # "selective" or "hemolytic"
    reference: str
    safety_score: float
    selectivity_proxy: float
    hydrophobic_fraction: float
    hydrophobic_moment: float
    charge_density_ph74: float
    gravy: float
    label: int  # 0 = selective, 1 = hemolytic (positive = dangerous)


@dataclass
class SelectivityResult:
    n_selective: int
    n_hemolytic: int
    safety_auroc: float
    selectivity_proxy_auroc: float
    hydrophobic_fraction_auroc: float
    hydrophobic_moment_auroc: float
    charge_density_auroc: float
    gravy_auroc: float
    rows: list[SelectivityRow] = field(default_factory=list)
    blind_spots: list[str] = field(default_factory=list)
    verdict: str = ""
    disclaimer: str = ""

    def to_dict(self) -> dict:
        return {
            "n_selective": self.n_selective,
            "n_hemolytic": self.n_hemolytic,
            "safety_auroc": self.safety_auroc,
            "selectivity_proxy_auroc": self.selectivity_proxy_auroc,
            "hydrophobic_fraction_auroc": self.hydrophobic_fraction_auroc,
            "hydrophobic_moment_auroc": self.hydrophobic_moment_auroc,
            "charge_density_auroc": self.charge_density_auroc,
            "gravy_auroc": self.gravy_auroc,
            "blind_spots": self.blind_spots,
            "verdict": self.verdict,
            "disclaimer": self.disclaimer,
            "per_peptide": [
                {
                    "id": r.id,
                    "sequence": r.sequence,
                    "family": r.family,
                    "hemolysis_class": r.hemolysis_class,
                    "safety_score": r.safety_score,
                    "selectivity_proxy": r.selectivity_proxy,
                    "hydrophobic_fraction": r.hydrophobic_fraction,
                    "hydrophobic_moment": r.hydrophobic_moment,
                    "charge_density_ph74": r.charge_density_ph74,
                    "gravy": r.gravy,
                }
                for r in self.rows
            ],
        }


def _auc_wilcoxon(pos_scores: list[float], neg_scores: list[float]) -> float:
    """AUROC via Wilcoxon-Mann-Whitney. Reuses the same method as retrospective.py.

    For selectivity: 'positive' = hemolytic (label=1, the dangerous class we want to detect).
    A good safety score should be LOW for hemolytic peptides, so AUROC is computed as
    P(safety_selective > safety_hemolytic) — i.e., safety_score should rank selectives higher.
    We achieve this by treating selective as 'pos' and hemolytic as 'neg' for safety_score.
    For other features (hydrophobic_fraction, gravy, hydrophobic_moment), higher values
    correlate with hemolysis, so we treat hemolytic as 'pos' and selective as 'neg'.
    """
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


def run_selectivity_benchmark(
    panel_csv: str | Path = "examples/validation/selectivity_panel.csv",
) -> SelectivityResult:
    """Run the selectivity benchmark on the reference panel.

    Returns a SelectivityResult with AUROC for each feature and per-peptide scores.
    """
    from openamp_foundry.features.physchem import compute_features
    from openamp_foundry.scoring.safety import safety_score

    rows: list[SelectivityRow] = []

    with open(panel_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            seq = row["sequence"].strip().upper()
            features = compute_features(seq)
            safe = safety_score(features)
            sel_proxy = features.get("selectivity_proxy", 1.0)

            hemolysis_class = row["hemolysis_class"].strip().lower()
            label = 1 if hemolysis_class == "hemolytic" else 0

            rows.append(SelectivityRow(
                id=row["id"],
                sequence=seq,
                family=row["family"],
                hemolysis_class=hemolysis_class,
                reference=row["reference"],
                safety_score=safe,
                selectivity_proxy=sel_proxy,
                hydrophobic_fraction=features["hydrophobic_fraction"],
                hydrophobic_moment=features.get("hydrophobic_moment", 0.0),
                charge_density_ph74=features.get("charge_density_ph74", 0.0),
                gravy=features.get("gravy", 0.0),
                label=label,
            ))

    selective = [r for r in rows if r.label == 0]
    hemolytic = [r for r in rows if r.label == 1]

    # For safety_score and selectivity_proxy: higher = safer = should rank selectives higher.
    # AUROC = P(score_selective > score_hemolytic).
    sel_safety = [r.safety_score for r in selective]
    hly_safety = [r.safety_score for r in hemolytic]
    safety_auroc = round(_auc_wilcoxon(sel_safety, hly_safety), 4)

    sel_proxy = [r.selectivity_proxy for r in selective]
    hly_proxy = [r.selectivity_proxy for r in hemolytic]
    proxy_auroc = round(_auc_wilcoxon(sel_proxy, hly_proxy), 4)

    # For hydrophobic_fraction, hydrophobic_moment, gravy: higher = more hemolytic.
    # AUROC = P(score_hemolytic > score_selective).
    hly_hfrac = [r.hydrophobic_fraction for r in hemolytic]
    sel_hfrac = [r.hydrophobic_fraction for r in selective]
    hfrac_auroc = round(_auc_wilcoxon(hly_hfrac, sel_hfrac), 4)

    hly_muh = [r.hydrophobic_moment for r in hemolytic]
    sel_muh = [r.hydrophobic_moment for r in selective]
    muh_auroc = round(_auc_wilcoxon(hly_muh, sel_muh), 4)

    hly_gravy = [r.gravy for r in hemolytic]
    sel_gravy = [r.gravy for r in selective]
    gravy_auroc = round(_auc_wilcoxon(hly_gravy, sel_gravy), 4)

    # Charge density: ambiguous — high charge can mean either strong AMP (selective)
    # or non-specific disruption. Test both directions and report the one that
    # discriminates (higher = hemolytic).
    hly_charge = [r.charge_density_ph74 for r in hemolytic]
    sel_charge = [r.charge_density_ph74 for r in selective]
    charge_auroc = round(_auc_wilcoxon(hly_charge, sel_charge), 4)

    # Identify blind spots: hemolytic peptides with safety_score >= 0.90
    blind_spots: list[str] = []
    for r in hemolytic:
        if r.safety_score >= 0.90:
            blind_spots.append(
                f"{r.id} ({r.family}): safety_score={r.safety_score} but "
                f"literature-hemolytic — known blind spot"
            )
    # Also flag selective peptides with low safety_score (false positives)
    for r in selective:
        if r.safety_score < 0.70:
            blind_spots.append(
                f"{r.id} ({r.family}): safety_score={r.safety_score} but "
                f"literature-selective — false alarm"
            )

    # Verdict: safety_score should achieve AUROC > 0.65 to be considered useful
    # for selectivity discrimination. Below 0.55 means it's worse than random.
    if safety_auroc >= 0.70:
        verdict = "STRONG — safety_score discriminates selective from hemolytic"
    elif safety_auroc >= 0.65:
        verdict = "MODERATE — safety_score has partial selectivity signal"
    elif safety_auroc >= 0.55:
        verdict = "WEAK — safety_score barely above random for selectivity"
    else:
        verdict = "FAILED — safety_score does NOT discriminate selective from hemolytic"

    disclaimer = (
        "This benchmark tests whether physicochemical proxies correlate with "
        "literature hemolysis classifications on a small reference panel (n=18). "
        "It does NOT predict hemolysis for novel sequences. Wet-lab hemolysis "
        "assay remains mandatory for all candidates. The panel is small and "
        "class-imbalanced; results are indicative, not definitive."
    )

    return SelectivityResult(
        n_selective=len(selective),
        n_hemolytic=len(hemolytic),
        safety_auroc=safety_auroc,
        selectivity_proxy_auroc=proxy_auroc,
        hydrophobic_fraction_auroc=hfrac_auroc,
        hydrophobic_moment_auroc=muh_auroc,
        charge_density_auroc=charge_auroc,
        gravy_auroc=gravy_auroc,
        rows=rows,
        blind_spots=blind_spots,
        verdict=verdict,
        disclaimer=disclaimer,
    )
