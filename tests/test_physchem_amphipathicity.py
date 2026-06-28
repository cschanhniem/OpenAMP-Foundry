"""Tests for hydrophobic moment (amphipathicity) feature."""
from __future__ import annotations

import pytest

from openamp_foundry.features.physchem import (
    compute_features,
    hydrophobic_moment,
    max_windowed_hydrophobic_moment,
)


class TestHydrophobicMoment:
    def test_empty_sequence_returns_zero(self):
        assert hydrophobic_moment("") == 0.0

    def test_uniform_sequence_low_moment(self):
        # All same amino acid → sine/cosine terms distribute evenly → low moment
        result = hydrophobic_moment("AAAAAAAAAA")
        assert isinstance(result, float)
        assert result >= 0.0

    def test_alternating_hydrophobic_polar_has_higher_moment(self):
        # Alternating hydrophobic/polar gives high periodicity
        # e.g. KLKLKLKL at 100deg/residue should show amphipathic character
        # vs KKKKKKKK which is all charged (low hydrophobic contrast)
        seq_amphipathic = "KLKLKLKL"
        seq_uniform = "KKKKKKKK"
        result_amph = hydrophobic_moment(seq_amphipathic)
        result_unif = hydrophobic_moment(seq_uniform)
        assert result_amph > result_unif, (
            f"Amphipathic KLKLKLKL moment ({result_amph:.4f}) should exceed "
            f"uniform KKKKKKKK moment ({result_unif:.4f})"
        )

    def test_known_amp_has_nonzero_moment(self):
        # KWKLFKKIGAVLKVL is a classic AMP (magainin analogue)
        result = hydrophobic_moment("KWKLFKKIGAVLKVL")
        assert result > 0.0

    def test_returns_float_rounded_to_4dp(self):
        result = hydrophobic_moment("KWKLFKK")
        assert isinstance(result, float)
        # Check 4 decimal places
        assert result == round(result, 4)

    def test_single_residue(self):
        result = hydrophobic_moment("K")
        # sin(0) = 0, cos(0) = 1 → moment = |H_K * cos(0)| / 1 = |H_K|
        assert isinstance(result, float)


class TestComputeFeaturesAmphipathicity:
    def test_hydrophobic_moment_in_features(self):
        features = compute_features("KWKLFKKIGAVLKVL")
        assert "hydrophobic_moment" in features
        assert isinstance(features["hydrophobic_moment"], float)
        assert features["hydrophobic_moment"] >= 0.0

    def test_empty_sequence_does_not_crash(self):
        features = compute_features("")
        assert "hydrophobic_moment" in features
        assert features["hydrophobic_moment"] == 0.0

    def test_all_canonical_amino_acids(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        features = compute_features(seq)
        assert "hydrophobic_moment" in features
        assert features["hydrophobic_moment"] >= 0.0


class TestMaxWindowedHydrophobicMoment:
    """Tests for the Eisenberg-standard windowed mu_h (window=11)."""

    def test_empty_sequence_returns_zero(self):
        assert max_windowed_hydrophobic_moment("") == 0.0

    def test_short_sequence_equals_full_seq_mu_h(self):
        # For sequences <= 11 AA, windowed == full-sequence (only one window)
        seq = "KWKLFKK"  # 7 residues < 11
        windowed = max_windowed_hydrophobic_moment(seq, window=11)
        full = hydrophobic_moment(seq)
        assert windowed == pytest.approx(full, abs=1e-3)

    def test_long_sequence_windowed_ge_full_seq(self):
        # For sequences > 11 AA, windowed >= full-sequence by definition
        # (best 11-residue window >= average over all residues)
        seq = "GIGKFLHSAKKFGKAFVGEIMNS"  # magainin-2, 23 residues
        windowed = max_windowed_hydrophobic_moment(seq, window=11)
        full = hydrophobic_moment(seq)
        assert windowed >= full - 1e-4, (
            f"Windowed mu_h ({windowed:.4f}) should be >= full-seq mu_h ({full:.4f})"
        )

    def test_magainin2_windowed_substantially_higher_than_full(self):
        # Magainin-2: windowed mu_h should be meaningfully higher than full-seq
        # (helical segment is concentrated; full-seq is diluted by terminal coil)
        magainin2 = "GIGKFLHSAKKFGKAFVGEIMNS"
        windowed = max_windowed_hydrophobic_moment(magainin2, window=11)
        full = hydrophobic_moment(magainin2)
        # Empirically: windowed ~0.68, full ~0.45 → ratio ~1.5×
        assert windowed > full * 1.2, (
            f"Magainin-2 windowed ({windowed:.4f}) should be >1.2× full-seq ({full:.4f})"
        )

    def test_returns_float_rounded_to_4dp(self):
        result = max_windowed_hydrophobic_moment("KWKLFKKIGAVLKVL")
        assert isinstance(result, float)
        assert result == round(result, 4)

    def test_single_residue_equals_full_moment(self):
        result = max_windowed_hydrophobic_moment("K", window=11)
        full = hydrophobic_moment("K")
        assert result == pytest.approx(full, abs=1e-3)

    def test_window_equal_to_sequence_length(self):
        seq = "KWKLFKKIGAV"  # exactly 11 residues
        windowed = max_windowed_hydrophobic_moment(seq, window=11)
        full = hydrophobic_moment(seq)
        assert windowed == pytest.approx(full, abs=1e-3)

    def test_nonnegative(self):
        for seq in ["", "A", "KWKLFKK", "GIGKFLHSAKKFGKAFVGEIMNS", "EEIEIEIEIEIEIEE"]:
            result = max_windowed_hydrophobic_moment(seq)
            assert result >= 0.0, f"max_windowed_hydrophobic_moment({seq!r}) returned {result}"

    def test_custom_window(self):
        seq = "KWKLFKKIGAVLKVL"
        w7 = max_windowed_hydrophobic_moment(seq, window=7)
        w11 = max_windowed_hydrophobic_moment(seq, window=11)
        # Both must be non-negative; neither is bounded at 1.0 (Eisenberg values can
        # be large for tightly packed hydrophobic windows — clamping happens in the scorer)
        assert w7 >= 0.0
        assert w11 >= 0.0


class TestComputeFeaturesMaxMuH:
    """Tests that compute_features() includes max_hydrophobic_moment."""

    def test_max_hydrophobic_moment_key_present(self):
        features = compute_features("KWKLFKKIGAVLKVL")
        assert "max_hydrophobic_moment" in features, (
            "compute_features() must return 'max_hydrophobic_moment'"
        )

    def test_max_hydrophobic_moment_is_float(self):
        features = compute_features("KWKLFKKIGAVLKVL")
        assert isinstance(features["max_hydrophobic_moment"], float)

    def test_max_hydrophobic_moment_nonnegative(self):
        for seq in ["", "A", "KWKLFKK", "GIGKFLHSAKKFGKAFVGEIMNS"]:
            features = compute_features(seq)
            assert features["max_hydrophobic_moment"] >= 0.0

    def test_max_hydrophobic_moment_ge_hydrophobic_moment_for_known_amps(self):
        # For known helical AMPs with a concentrated amphipathic segment, the windowed
        # value should be >= the full-seq value (good amphipathic signal).
        # NOTE: this invariant does NOT hold for non-amphipathic or uniform sequences
        # (see test_windowed_can_be_less_than_full_for_uniform_sequences below).
        for seq in ["KWKLFKKIGAVLKVL", "GIGKFLHSAKKFGKAFVGEIMNS", "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"]:
            features = compute_features(seq)
            assert features["max_hydrophobic_moment"] >= features["hydrophobic_moment"] - 1e-4, (
                f"max_mu_h ({features['max_hydrophobic_moment']:.4f}) < mu_h "
                f"({features['hydrophobic_moment']:.4f}) for {seq}"
            )

    def test_windowed_can_be_less_than_full_for_uniform_sequences(self):
        # For non-amphipathic sequences (e.g. all-Ala, all-Glu), the windowed mu_h
        # (normalised by window=11) can be LESS than full-seq mu_h (normalised by len>11).
        # activity_likeness_score() handles this correctly via max(); direct callers
        # must not assume max_hydrophobic_moment >= hydrophobic_moment.
        seq = "AAAAAAAAAAAA"  # 12-A, uniform → no concentrated amphipathic segment
        windowed = max_windowed_hydrophobic_moment(seq, window=11)
        full = hydrophobic_moment(seq)
        assert windowed < full, (
            f"For uniform sequence {seq!r}, expected windowed ({windowed:.4f}) < "
            f"full-seq ({full:.4f}) to document the non-invariant"
        )

    def test_short_seq_max_equals_full(self):
        features = compute_features("KWKLFKK")  # 7 residues, window=11 → only one window
        assert features["max_hydrophobic_moment"] == pytest.approx(
            features["hydrophobic_moment"], abs=1e-3
        )
