from __future__ import annotations


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def serum_stability_score(features: dict) -> float:
    """Predict relative serum stability from interior protease-site density.

    Short cationic AMPs are degraded primarily by serum trypsin and chymotrypsin.
    Each interior K/R site adds trypsin susceptibility; F/W/Y sites add chymotrypsin
    susceptibility.  Higher site density → lower predicted stability.

    Literature basis:
    - Hilpert et al. (2006) J Antimicrob Chemother: short cationic peptides with
      ≥3 interior K/R sites typically have serum t½ < 30 min.
    - Wade et al. (1990) PNAS: D-amino acid substitution extends t½ by >10×.

    Returns a value in [0, 1]:
      1.0 = no interior cleavage sites (e.g., pure alanine or hydrophobic peptide)
      0.5 = moderate site density (~0.25 per residue)
      0.0 = very high site density (≥0.5 per residue)

    This score is informational — it does not currently gate the ensemble.
    Use it to stratify candidates by predicted serum longevity and to flag
    candidates that require stability engineering before clinical translation.
    """
    length = features.get("length", 0)
    if length == 0:
        return 1.0

    trypsin_density = features.get("trypsin_site_density", 0.0)
    n_chymo = features.get("interior_chymotrypsin_sites", 0)
    chymo_density = round(n_chymo / length, 4)

    # Trypsin is the primary serum protease for cationic peptides; weight 2:1 vs chymotrypsin
    combined_density = (2.0 * trypsin_density + 1.0 * chymo_density) / 3.0

    # At combined_density = 0.5 (every other residue is a cleavage site), score → 0
    score = 1.0 - clamp01(combined_density / 0.5)
    return round(score, 4)
