"""Tests for retrospective AUROC benchmark module."""
from __future__ import annotations

import csv

import pytest

from openamp_foundry.benchmark.retrospective import (
    _auc_wilcoxon,
    _recall_at_k,
    run_retrospective_benchmark,
)


class TestAucWilcoxon:
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

    def test_empty_pos_returns_half(self):
        assert _auc_wilcoxon([], [0.5, 0.6]) == pytest.approx(0.5)

    def test_ties_counted_as_half(self):
        pos = [0.6]
        neg = [0.6]
        assert _auc_wilcoxon(pos, neg) == pytest.approx(0.5)


class TestRecallAtK:
    def test_all_positives_at_top(self):
        labels = [1, 1, 1, 0, 0]
        assert _recall_at_k(labels, k=3) == pytest.approx(1.0)

    def test_no_positives_at_top(self):
        labels = [0, 0, 0, 1, 1]
        assert _recall_at_k(labels, k=3) == pytest.approx(0.0)

    def test_half_recall(self):
        labels = [1, 0, 1, 0, 0, 0]
        # 1 of 2 positives in top-2
        assert _recall_at_k(labels, k=2) == pytest.approx(0.5)

    def test_no_positives_returns_zero(self):
        assert _recall_at_k([0, 0, 0], k=2) == pytest.approx(0.0)


class TestRunRetrospectiveBenchmark:
    @pytest.fixture
    def mini_amp_csv(self, tmp_path):
        p = tmp_path / "amps.csv"
        rows = [
            {"id": "AMP-001", "sequence": "KWKLFKKIGAVLKVL", "family": "template", "reference": "test", "label": 1},
            {"id": "AMP-002", "sequence": "RRWQWRMKKLG", "family": "rrw", "reference": "test", "label": 1},
            {"id": "AMP-003", "sequence": "GIGKFLHSAKKFGKAFVGEIMNS", "family": "magainin", "reference": "test", "label": 1},
        ]
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "sequence", "family", "reference", "label"])
            w.writeheader()
            w.writerows(rows)
        return p

    @pytest.fixture
    def mini_decoy_csv(self, tmp_path):
        p = tmp_path / "decoys.csv"
        # These are composition-shuffled (all-G/all-A decoys that score low)
        rows = [
            {"id": "DECOY-001", "sequence": "GGGGGGGGGGGGGGG", "family": "shuffled", "source_id": "AMP-001", "label": 0},
            {"id": "DECOY-002", "sequence": "AAAAAAAAAAA", "family": "shuffled", "source_id": "AMP-002", "label": 0},
            {"id": "DECOY-003", "sequence": "GGGGGGGGGGGGGGGGGGGGGGG", "family": "shuffled", "source_id": "AMP-003", "label": 0},
        ]
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "sequence", "family", "source_id", "label"])
            w.writeheader()
            w.writerows(rows)
        return p

    def test_returns_required_keys(self, mini_amp_csv, mini_decoy_csv):
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        for key in ["auroc", "n_positives", "n_negatives", "interpretation", "top_ranked", "disclaimer"]:
            assert key in result, f"Missing key: {key}"

    def test_auroc_between_0_and_1(self, mini_amp_csv, mini_decoy_csv):
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        assert 0.0 <= result["auroc"] <= 1.0

    def test_amps_outrank_all_g_decoys(self, mini_amp_csv, mini_decoy_csv):
        # Known AMPs should clearly outrank all-G / all-A decoys
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        assert result["auroc"] > 0.7, (
            f"Known AMPs should strongly outrank degenerate decoys (AUROC={result['auroc']:.4f})"
        )

    def test_n_counts_correct(self, mini_amp_csv, mini_decoy_csv):
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        assert result["n_positives"] == 3
        assert result["n_negatives"] == 3
        assert result["n_total"] == 6

    def test_recall_keys_present(self, mini_amp_csv, mini_decoy_csv):
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        assert "recall_at_10" in result or "recall_at_3" in result or any(
            k.startswith("recall_at") for k in result
        )

    def test_interpretation_is_string(self, mini_amp_csv, mini_decoy_csv):
        result = run_retrospective_benchmark(mini_amp_csv, mini_decoy_csv)
        assert isinstance(result["interpretation"], str)
        assert len(result["interpretation"]) > 10

    def test_full_benchmark_with_real_data(self):
        """Run on the actual validation dataset and report the AUROC honestly."""
        from pathlib import Path
        amp_csv = Path("examples/validation/known_amps.csv")
        decoy_csv = Path("examples/validation/scrambled_decoys.csv")
        if not amp_csv.exists() or not decoy_csv.exists():
            pytest.skip("Validation data not found — run from project root")
        result = run_retrospective_benchmark(amp_csv, decoy_csv)
        # We do NOT assert a specific AUROC — we report what the model actually achieves.
        # The benchmark result determines whether to proceed with wet-lab synthesis.
        assert 0.0 <= result["auroc"] <= 1.0
        # Known AMPs must outrank at least random (> 0.50) — otherwise the model is broken
        assert result["auroc"] > 0.50, (
            f"AUROC={result['auroc']:.4f}: model performs below random on known AMPs vs shuffled decoys. "
            "The scoring model is broken and must not be used for candidate nomination."
        )
        # Print for visibility
        print(f"\nRetrospective AUROC: {result['auroc']:.4f}")
        print(f"Interpretation: {result['interpretation']}")
        print(f"Recall@10: {result.get('recall_at_10', 'N/A')}")
        print(f"Recall@20: {result.get('recall_at_20', 'N/A')}")
