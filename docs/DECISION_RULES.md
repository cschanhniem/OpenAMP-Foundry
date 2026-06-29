# Pipeline Decision Rules

> **Purpose:** Pre-registered pass/fail gates for each stage of the OpenAMP Foundry
> pipeline. Thresholds are hardcoded (see `src/openamp_foundry/gates/gate_checker.py`)
> to prevent cherry-picking after inspecting results.
>
> **Locked:** 2026-06-29 (Sprint 5)

---

## Gate 1 — AUROC Benchmark

| Field | Value |
|-------|-------|
| Threshold | AUROC ≥ **0.70** |
| Measure | Pipeline AUROC on expanded 95+96 benchmark |
| Config | `configs/pipeline.yaml` |
| Action | If FAIL: Do not proceed to synthesis |

## Gate 2 — Leakage Guard

| Field | Value |
|-------|-------|
| Threshold | Recall@43 < **0.60** |
| Measure | Fraction of positives recovered in top 43 ranked |
| Rationale | Recall@43 > 0.60 suggests near-duplicate memorisation |
| Action | If FAIL: Investigate reference-set contamination |

## Gate 3 — Model Disagreement

| Field | Value |
|-------|-------|
| Threshold | Top-3 candidate |activity − boman_activity| < **0.45** |
| Measure | Absolute difference between two independent scorers |
| Config | `configs/pipeline.yaml`, `configs/phase3.yaml` |
| Action | If FAIL: Review scorer weights or threshold |

## Gate 4 — Top-10 Recall

| Field | Value |
|-------|-------|
| Threshold | Recall@10 > **0.0** |
| Measure | At least 1 positive in top 10 ranked candidates |
| Action | If FAIL: Scoring model likely broken |

## Gate 5 — Interpretation

| Field | Value |
|-------|-------|
| Threshold | Benchmark interpretation must be **STRONG** |
| Measure | String match of `interpretation` field |
| Action | If FAIL: Do not proceed |

## Gate 6 — External Predictor Consensus

| Field | Value |
|-------|-------|
| Threshold | ≥3/5 external tools positive |
| Measure | Manual web submission results |
| Status | PENDING (see `outputs/external_predict_checklist.md`) |
| Action | If FAIL: Expert reviewer override required |

## Gate 7 — Human Expert Review

| Field | Value |
|-------|-------|
| Threshold | Expert reviewer APPROVE or CONDITIONAL |
| Measure | Completed reviewer questionnaire |
| Status | PENDING (see `outputs/questionnaire/`) |
| Action | If REJECT: Exclude from synthesis |

---

## CLI Usage

```bash
# Check all gates:
openamp-foundry gate-check \
    --validation-json outputs/validate_scoring_report.json

# Check specific gate:
openamp-foundry gate-check --gate 1 \
    --validation-json outputs/validate_scoring_report.json
```

## Hardcoding Rule

All gate thresholds in `gates/gate_checker.py` are **hardcoded** (not configurable
via YAML). To change a threshold: modify source code, commit, and create a new release.
This prevents silent threshold drift that could mask performance regression.
