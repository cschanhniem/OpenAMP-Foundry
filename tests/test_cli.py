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
    assert "MET" in text or "oxidation" in text.lower()


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
