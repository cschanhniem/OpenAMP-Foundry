# External Predictor Consensus

## Purpose

Aggregate predictions from 5 independent AMP prediction tools into a single
consensus verdict per candidate. A candidate with ≥3/5 tool agreement is
classified as **CONFIDENT** and recommended for synthesis.

## Tool Summary

| Tool | Method | Positive label | URL |
|------|--------|---------------|-----|
| CAMPR4 | SVM + RF + ANN + DT ensemble | AMP | camp3.bicnirrh.res.in |
| AMPScanner v2 | LSTM | Antimicrobial | dveltri.com/ascan/v2 |
| dbAMP 2.0 | Random Forest | AMP | awi.cuhk.edu.cn/dbAMP |
| AntiCP 2.0 | SVM | ACP (not AMP) | webs.iiitd.edu.in/raghava/anticp2 |
| Macrel | Random Forest | AMP | big-data-biology.org/software/macrel |

## Consensus Rules

| Agreement | Verdict | Action |
|:---------:|---------|--------|
| ≥3/5 | CONFIDENT | Proceed to synthesis |
| 2/5 | UNCERTAIN | Expert review required |
| ≤1/5 | WEAK | Do not synthesise |

## Current Status

**PENDING** — Results not yet obtained. All 5 tools require manual web submission.
See `outputs/external_predict_checklist.md` for the submission guide.

## CLI Usage

```bash
# After filling in results CSV:
make external-consensus RESULTS=outputs/external_predict_results.csv

# Or manually:
openamp-foundry external-consensus \
    --pilot-csv outputs/pilot_panel.csv \
    --results-csv outputs/external_predict_results.csv \
    --out outputs/external_consensus_report.md
```

## See Also

- `outputs/external_predict_checklist.md` — Manual submission guide
- `src/openamp_foundry/reports/external_consensus.py` — Implementation
- `docs/ASSAY_PREREGISTRATION.md` — Pre-registered assay protocol
