"""CLI integration tests."""
from __future__ import annotations

import json

from openamp_foundry.cli import main


def test_rank_command_success(tmp_path):
    out = str(tmp_path / "ranked.jsonl")
    ret = main([
        "rank",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--out", out,
    ])
    assert ret == 0


def test_rank_command_with_report_and_certs(tmp_path):
    out = str(tmp_path / "ranked.jsonl")
    report = str(tmp_path / "report.md")
    certs = str(tmp_path / "certs")
    ret = main([
        "rank",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--out", out,
        "--report", report,
        "--cert-dir", certs,
    ])
    assert ret == 0


def test_validate_command_success(tmp_path):
    # First generate a certificate
    out = str(tmp_path / "ranked.jsonl")
    certs = str(tmp_path / "certs")
    main([
        "rank",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--out", out,
        "--cert-dir", certs,
    ])
    cert_files = list((tmp_path / "certs").glob("*.json"))
    assert cert_files, "No certificates were generated"
    ret = main([
        "validate",
        "--certificate", str(cert_files[0]),
        "--schema", "schemas/candidate.schema.json",
    ])
    assert ret == 0


def test_bench_leakage_detects_duplicates(tmp_path, capsys):
    # Demo candidates 1, 2, 5 are exact copies of references
    ret = main([
        "bench", "leakage",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
    ])
    assert ret == 0
    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["near_duplicate_count"] == 3
    assert result["warning"] is not None


def test_bench_leakage_no_duplicates(tmp_path, capsys):
    # Use negative examples as candidates — they won't match the reference AMPs
    ret = main([
        "bench", "leakage",
        "--candidates", "examples/negative/demo_negative_peptides.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--threshold", "0.90",
    ])
    assert ret == 0
    captured = capsys.readouterr()
    result = json.loads(captured.out)
    assert result["near_duplicate_count"] == 0
    assert result["warning"] is None


def test_bench_leakage_output_file(tmp_path, capsys):
    out = str(tmp_path / "leakage_report.json")
    ret = main([
        "bench", "leakage",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--out", out,
    ])
    assert ret == 0
    data = json.loads((tmp_path / "leakage_report.json").read_text())
    assert "near_duplicates" in data


def test_report_contains_disclaimer(tmp_path):
    out = str(tmp_path / "ranked.jsonl")
    report = str(tmp_path / "report.md")
    main([
        "rank",
        "--candidates", "examples/sequences/demo_candidates.csv",
        "--references", "examples/known_reference/demo_known_amps.csv",
        "--out", out,
        "--report", report,
    ])
    text = (tmp_path / "report.md").read_text()
    assert "NOT validated biological predictors" in text
    assert "no antimicrobial activity has been demonstrated" in text.lower() or "No antimicrobial activity" in text


def test_presynth_qc_command_returns_zero(tmp_path):
    panel = tmp_path / "panel.csv"
    panel.write_text(
        "candidate_id,sequence,source\n"
        "SEED-001,KWKLFKKIGAVLKVL,test\n"
        "SEED-002,RRWQWRMKKLG,test\n"
    )
    out = str(tmp_path / "report.md")
    ret = main(["presynth-qc", "--panel-csv", str(panel), "--out", out])
    assert ret == 0


def test_presynth_qc_command_creates_report(tmp_path):
    panel = tmp_path / "panel.csv"
    panel.write_text(
        "candidate_id,sequence,source\n"
        "SEED-001,KWKLFKKIGAVLKVL,test\n"
    )
    out_path = tmp_path / "qc_report.md"
    main(["presynth-qc", "--panel-csv", str(panel), "--out", str(out_path)])
    assert out_path.exists()
    text = out_path.read_text()
    assert "Pre-Synthesis QC Report" in text
    assert "SEED-001" in text


def test_presynth_qc_command_flags_met_residue(tmp_path):
    panel = tmp_path / "panel.csv"
    panel.write_text(
        "candidate_id,sequence,source\n"
        "MET-001,KRLMKKIGSAIKFL,test\n"
    )
    out_path = tmp_path / "qc_report.md"
    main(["presynth-qc", "--panel-csv", str(panel), "--out", str(out_path)])
    text = out_path.read_text()
    # Both the candidate ID and the MET flag must appear — avoids a false pass from
    # a generic preamble sentence that happens to contain the word "oxidation".
    assert "MET-001" in text and "MET" in text


def test_presynth_qc_command_contains_summary_table(tmp_path):
    panel = tmp_path / "panel.csv"
    panel.write_text(
        "candidate_id,sequence,source\n"
        "A,KWKLFKKIGAVLKVL,test\n"
        "B,AAAAAAAAGGGGGGGG,test\n"
    )
    out_path = tmp_path / "qc_report.md"
    main(["presynth-qc", "--panel-csv", str(panel), "--out", str(out_path)])
    text = out_path.read_text()
    assert "Summary Table" in text
    assert "Candidates checked: 2" in text


def test_validate_scoring_stdout_includes_n_and_auprc(tmp_path, capsys):
    """validate-scoring stdout must include n_positives, n_negatives, benchmark_type, auprc."""
    import pathlib
    amp_csv = pathlib.Path("examples/validation/known_amps.csv")
    bg_csv = pathlib.Path("examples/validation/random_background.csv")
    if not amp_csv.exists() or not bg_csv.exists():
        import pytest
        pytest.skip("Validation data not found — run from project root")
    out = str(tmp_path / "report.json")
    main([
        "validate-scoring",
        "--amp-csv", str(amp_csv),
        "--decoy-csv", str(bg_csv),
        "--config", "configs/pipeline.yaml",
        "--benchmark-type", "standard",
        "--out", out,
    ])
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert "n_positives" in data, "stdout missing n_positives"
    assert "n_negatives" in data, "stdout missing n_negatives"
    assert "benchmark_type" in data, "stdout missing benchmark_type"
    assert "auprc" in data, "stdout missing auprc"
    assert data["n_positives"] == 43
    assert data["n_negatives"] == 44
    assert data["benchmark_type"] == "standard"
    assert 0.0 < data["auprc"] < 1.0


def test_pilot_panel_malformed_jsonl_returns_error(tmp_path, capsys):
    """pilot-panel must return structured error on malformed JSONL, not crash."""
    bad_jsonl = tmp_path / "bad.jsonl"
    bad_jsonl.write_text('{"candidate_id": "X", "selected": true}\n{BROKEN JSON LINE\n')
    out_csv = str(tmp_path / "panel.csv")
    out_md = str(tmp_path / "panel.md")
    rc = main([
        "pilot-panel",
        "--ranked", str(bad_jsonl),
        "--out-csv", out_csv,
        "--out-md", out_md,
    ])
    assert rc == 1, "Expected non-zero exit code on malformed JSONL"
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["status"] == "error"
    assert "line 2" in data["message"]
    assert "line_preview" in data


def test_synthesis_order_missing_columns_returns_error(tmp_path, capsys):
    """synthesis-order must return structured error when panel CSV lacks required columns."""
    bad_panel = tmp_path / "bad_panel.csv"
    bad_panel.write_text("pilot_rank,candidate_id\n1,SEED-003_VAR_001\n")  # missing 'sequence'
    out_csv = str(tmp_path / "order.csv")
    rc = main([
        "synthesis-order",
        "--panel-csv", str(bad_panel),
        "--out-csv", out_csv,
    ])
    assert rc == 1, "Expected non-zero exit code on missing columns"
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["status"] == "error"
    assert "sequence" in data["message"]
