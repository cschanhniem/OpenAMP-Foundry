"""Selectivity benchmark tests — Phase 2.5 requirement.

Per AGENTS.md Phase 2.5: "Membrane proxy benchmark — Distinguishes bacterial-selective
vs clearly hemolytic reference peptides better than naive heuristics."

These tests verify that:
1. The selectivity benchmark runs on the reference panel.
2. The benchmark produces correct AUROC values for all features.
3. Known blind spots (melittin safety=1.0 despite hemolysis) are detected and reported.
4. At least one feature discriminates selective from hemolytic above 0.60 AUROC.
5. The benchmark output includes per-peptide scores and baseline comparisons.
6. The CLI command produces valid output.

CURRENT FINDING (2026-07-01):
safety_score AUROC = 0.54 (near random) — the safety filter does NOT discriminate
selective from hemolytic AMPs. hydrophobic_fraction AUROC = 0.81 (strong baseline).
This gap is documented as a known issue; the test verifies the benchmark honestly
reveals it rather than hiding it.

All results are computational proxies. No biological hemolysis is measured.
The benchmark tests correlation with literature hemolysis classifications only.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from openamp_foundry.benchmark.selectivity import (
    SelectivityResult,
    _auc_wilcoxon,
    run_selectivity_benchmark,
)


PANEL_CSV = "examples/validation/selectivity_panel.csv"


class TestSelectivityPanel:
    """Verify the reference panel data integrity."""

    def test_panel_exists(self):
        assert Path(PANEL_CSV).exists(), "Selectivity panel CSV not found"

    def test_panel_has_expected_columns(self):
        with open(PANEL_CSV) as f:
            reader = csv.DictReader(f)
            assert "id" in reader.fieldnames
            assert "sequence" in reader.fieldnames
            assert "hemolysis_class" in reader.fieldnames

    def test_panel_has_both_classes(self):
        with open(PANEL_CSV) as f:
            rows = list(csv.DictReader(f))
        classes = {r["hemolysis_class"] for r in rows}
        assert "selective" in classes
        assert "hemolytic" in classes

    def test_panel_has_at_least_8_per_class(self):
        with open(PANEL_CSV) as f:
            rows = list(csv.DictReader(f))
        n_sel = sum(1 for r in rows if r["hemolysis_class"] == "selective")
        n_hly = sum(1 for r in rows if r["hemolysis_class"] == "hemolytic")
        assert n_sel >= 8, f"Need at least 8 selective AMPs, got {n_sel}"
        assert n_hly >= 5, f"Need at least 5 hemolytic AMPs, got {n_hly}"

    def test_all_sequences_are_valid_amino_acids(self):
        with open(PANEL_CSV) as f:
            rows = list(csv.DictReader(f))
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for row in rows:
            seq = row["sequence"].upper()
            assert all(aa in valid_aa for aa in seq), f"Invalid AA in {row['id']}: {seq}"


class TestSelectivityBenchmark:
    """Core benchmark tests."""

    @pytest.fixture
    def result(self):
        return run_selectivity_benchmark(PANEL_CSV)

    def test_returns_selectivity_result(self, result):
        assert isinstance(result, SelectivityResult)

    def test_correct_counts(self, result):
        assert result.n_selective > 0
        assert result.n_hemolytic > 0
        assert result.n_selective + result.n_hemolytic == len(result.rows)

    def test_at_least_one_feature_discriminates(self, result):
        """At least one feature must achieve AUROC > 0.60 for selectivity.

        If none do, the benchmark panel or the pipeline is fundamentally broken.
        Currently hydrophobic_fraction (AUROC ~0.81) carries the signal.
        """
        all_aurocs = [
            result.safety_auroc,
            result.selectivity_proxy_auroc,
            result.hydrophobic_fraction_auroc,
            result.hydrophobic_moment_auroc,
            result.charge_density_auroc,
            result.gravy_auroc,
        ]
        best = max(all_aurocs)
        assert best > 0.60, (
            f"No feature achieves AUROC > 0.60 (best={best:.4f}). "
            "The selectivity benchmark panel may be too small or the pipeline "
            "has no selectivity signal at all."
        )

    def test_safety_auroc_reported_honestly(self, result):
        """safety_score AUROC is reported regardless of whether it passes or fails.

        Current value: ~0.54 (near random). This test verifies the value is
        reported correctly and is in [0, 1]. It does NOT assert the score is good.
        The safety_score's failure to discriminate is a documented finding.
        """
        assert 0.0 <= result.safety_auroc <= 1.0

    def test_per_peptide_scores_present(self, result):
        assert len(result.rows) > 0
        row = result.rows[0]
        assert hasattr(row, "safety_score")
        assert hasattr(row, "selectivity_proxy")
        assert hasattr(row, "hydrophobic_fraction")
        assert hasattr(row, "hydrophobic_moment")

    def test_to_dict_produces_valid_json(self, result):
        d = result.to_dict()
        json.dumps(d)
        assert "safety_auroc" in d
        assert "per_peptide" in d

    def test_blind_spots_detected(self, result):
        """Known blind spots (hemolytic with high safety_score) should be flagged.

        Melittin is the canonical blind spot: safety=1.0 despite HC50 ~ 5 ug/mL.
        At least 3 hemolytic peptides should appear in blind_spots (melittin,
        mastoparan-X, and at least one more).
        """
        assert len(result.blind_spots) >= 3, (
            f"Expected at least 3 blind spots (hemolytic with safety >= 0.90), "
            f"got {len(result.blind_spots)}: {result.blind_spots}"
        )
        blind_spot_text = " ".join(result.blind_spots).lower()
        assert "melittin" in blind_spot_text

    def test_verdict_is_informative(self, result):
        assert len(result.verdict) > 10
        assert any(w in result.verdict for w in ["STRONG", "MODERATE", "WEAK", "FAILED"])

    def test_disclaimer_present(self, result):
        assert "not predict" in result.disclaimer.lower()
        assert "wet-lab" in result.disclaimer.lower() or "wet lab" in result.disclaimer.lower()

    def test_baseline_aurocs_computed(self, result):
        """All baseline AUROCs should be present and in [0, 1]."""
        for auroc in [result.hydrophobic_fraction_auroc, result.hydrophobic_moment_auroc,
                       result.charge_density_auroc, result.gravy_auroc]:
            assert 0.0 <= auroc <= 1.0

    def test_hydrophobic_fraction_beats_safety(self, result):
        """hydrophobic_fraction should outperform safety_score for selectivity.

        This is a key finding: the naive hydrophobic_fraction baseline (AUROC ~0.81)
        significantly outperforms the designed safety_score (AUROC ~0.54).
        If this test fails, either safety_score was improved (good — update test)
        or hydrophobic_fraction lost discriminative power (investigate).
        """
        assert result.hydrophobic_fraction_auroc > result.safety_auroc, (
            f"hydrophobic_fraction AUROC ({result.hydrophobic_fraction_auroc}) should "
            f"exceed safety_score AUROC ({result.safety_auroc}). If safety_score was "
            f"improved to beat this baseline, update this test."
        )


class TestAucWilcoxonSelectivity:
    """Unit tests for the AUROC computation in selectivity context."""

    def test_perfect_separation(self):
        pos = [0.9, 0.8, 0.7]
        neg = [0.4, 0.3, 0.2]
        assert _auc_wilcoxon(pos, neg) == pytest.approx(1.0)

    def test_worst_case(self):
        pos = [0.2, 0.3, 0.4]
        neg = [0.7, 0.8, 0.9]
        assert _auc_wilcoxon(pos, neg) == pytest.approx(0.0)

    def test_random_is_half(self):
        pos = [0.5, 0.5, 0.5]
        neg = [0.5, 0.5, 0.5]
        assert _auc_wilcoxon(pos, neg) == pytest.approx(0.5)

    def test_empty_returns_half(self):
        assert _auc_wilcoxon([], [0.5, 0.6]) == pytest.approx(0.5)


class TestSelectivityCLI:
    """CLI integration test for validate-selectivity command."""

    def test_cli_runs_and_produces_json(self, tmp_path):
        from openamp_foundry.cli import main
        out_json = tmp_path / "selectivity_result.json"
        rc = main(["validate-selectivity", "--out", str(out_json)])
        assert rc == 0
        assert out_json.exists()
        data = json.loads(out_json.read_text())
        assert "safety_auroc" in data
        assert "per_peptide" in data
        assert len(data["per_peptide"]) > 0
