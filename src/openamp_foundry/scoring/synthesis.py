from __future__ import annotations

from openamp_foundry.scoring.activity import clamp01


def synthesis_feasibility_score(features: dict, valid_sequence: bool = True) -> float:
    """Estimate solid-phase synthesis difficulty [0, 1]; higher = easier to synthesise.

    Penalty sources:
    - Length > 30: longer chains accumulate deletion errors and solubility issues.
    - Length < 8: too short to be reliably purified / characterised.
    - Longest single-AA repeat run ≥ 5: coupling efficiency drops on homo-repeat stretches.
    - Cysteine fraction > 20%: disulphide scrambling and side-chain protection cost.
    - Aggregation propensity > 0: interior hydrophobic runs (VILMFW ≥ 4) and high
      beta-branched density (Val/Ile/Thr) cause on-resin aggregation and poor solubility.
      References: Quittot et al. (2017) Protein Sci; Wurth et al. (2006) J Mol Biol.
    - Proline fraction > 15%: N-methylated backbone requires extended activation (slow couplings
      at XP junctions); DKP formation risk at N-terminal dipeptides containing Pro during
      piperidine Fmoc deprotection; overall coupling efficiency reduced. References: Barlos et
      al. (1989) Int J Peptide Protein Res; Quibell et al. (1994) J Am Chem Soc; Fischer (2003)
      Curr Opin Drug Discov Devel.
    """
    if not valid_sequence:
        return 0.0
    length = features["length"]
    repeat_run = features["longest_repeat_run"]
    cys = features["cysteine_fraction"]
    agg = features.get("aggregation_propensity", 0.0)
    pro = features.get("proline_fraction", 0.0)

    score = 1.0
    if length > 30:
        score -= min((length - 30) * 0.04, 0.40)
    if length < 8:
        score -= 0.30
    if repeat_run >= 5:
        score -= 0.10
    if cys > 0.20:
        score -= 0.15
    if agg > 0.0:
        score -= min(agg * 0.25, 0.20)
    if pro > 0.15:
        score -= 0.10
    return round(clamp01(score), 4)
