# Benchmark Card — OpenAMP Foundry v0.5.x

> **Purpose:** Single-page summary of benchmark methodology, data, and metrics.
> **Last updated:** 2026-07-01 (selectivity benchmark added)

---

## Identity

| Field | Value |
|-------|-------|
| Pipeline version | v0.5.x |
| Benchmark type | Retrospective AUROC (composition-matched decoys) |
| Date | 2026-06-29 |
| Config | `configs/pipeline.yaml` (primary), `configs/phase3.yaml` (synthesis gate) |
| Commit | See `docs/METRICS_CURRENT.md` |

## Data

| Set | Count | Description |
|-----|:-----:|-------------|
| Positives | **95** | Public-domain AMPs from 12 taxonomic classes |
| Negatives (standard) | **96** | Length-matched random decoys (Swiss-Prot residue frequencies) |
| Negatives (strict) | **95** | Per-sequence composition-matched shuffles |
| Total (standard) | **191** | Expanded from original 87 (PR #110) |
| Reference library | **120** | Unified AMP library (deduplicated curated_72 + expanded_95) |

## Metrics

| Metric | Standard (pipeline) | Phase3 gate | Strict |
|--------|:-------------------:|:-----------:|:------:|
| **AUROC** | **0.7832** | **0.7448** | 0.4335 |
| **CI₉₅** | 0.72–0.84 | 0.68–0.81 | — |
| **AUPRC** | 0.8164 | 0.7933 | — |
| Baseline AUPRC | 0.4974 | 0.4974 | — |
| Recall@10 | 0.1053 | 0.1053 | — |
| Recall@20 | 0.2105 | 0.2105 | — |
| Recall@43 | 0.4211 | 0.4000 | — |
| Bootstrap | 2000 resamples | 2000 resamples | — |
| **Interpretation** | **STRONG** | **STRONG** | Below random (expected) |

## Method

AUROC computed via Wilcoxon-Mann-Whitney statistic (concordant-pair enumeration).
AUPRC via trapezoidal integration of precision-recall curve (pessimistic tie-breaking).
Confidence intervals: percentile bootstrap (2000 resamples, seed=0).

## Known Biases

| Bias | Impact | Mitigation |
|------|--------|------------|
| Helical-centric scorer | β-sheet AMPs (defensins) score below panel | Noted in METRICS_CURRENT.md |
| Melittin safety blind spot | Safety=1.0 despite strong hemolysis; 5/8 hemolytic AMPs score safety=1.0 (selectivity AUROC=0.54) | Hemolysis assay mandatory for all; safety_score redesign identified as improvement target |
| Composition-matched scrambled AUROC < 0.5 | Model relies on composition > order; fails strict order-sensitivity test | OpenAMP is an **evidence-ranking tool**, not a validated sequence-order activity predictor. Composition-scrambled test confirms the model captures AMP-like composition, not sequence-order features. |
| Near-seed generation only | Novel sequence space not explored | Acknowledged limitation in all docs |

## Classification

OpenAMP performs well on standard decoys (AUROC 0.7832) but **fails the strict composition-scrambled benchmark (AUROC 0.4335, below random)**. This means the model primarily detects AMP-like amino acid composition rather than sequence-order-dependent antimicrobial features. This is expected for a physicochemical heuristic scorer and is appropriate for its intended use: **triage and ranking**, not deep biological prediction.

## Selectivity Benchmark

> Tests whether safety features distinguish selective from hemolytic AMPs.
> Panel: 10 selective + 8 hemolytic AMPs (literature hemolysis classifications).
> Run: `make validate-selectivity`

| Feature | AUROC | Verdict |
|---------|:-----:|---------|
| safety_score | 0.54 | **FAILED** — near random |
| selectivity_proxy | 0.61 | Weak |
| hydrophobic_fraction | 0.81 | Strong naive baseline |
| GRAVY | 0.76 | Strong naive baseline |
| hydrophobic_moment | 0.61 | Weak |

**Finding:** safety_score does NOT discriminate selective from hemolytic AMPs.
5/8 hemolytic AMPs score safety=1.0 (melittin, mastoparan-X, PMAP-23, bombolitin-II, polybia-MP1).
hydrophobic_fraction is the best single discriminator. Improvement target identified.

## Historical Baseline

| Point | Set | Pipeline AUROC | Phase3 AUROC |
|-------|:---:|:--------------:|:------------:|
| Current (PR #110) | 95+96 (n=191) | **0.7832** | 0.7448 |
| Pre-expansion (PR #72) | 43+44 (n=87) | 0.8420 | 0.8266 |
| Pre-face-bonus (PR #70) | 43+44 | 0.8348 | 0.8126 |
| Pre-windowed-mu_h (PR #66) | 43+44 | 0.8047 | 0.7846 |

## Change Log

| Date | Change | Author |
|------|--------|--------|
| 2026-07-01 | Selectivity benchmark added — safety_score AUROC=0.54 | OpenAMP CI |
| 2026-06-29 | Initial card — expanded benchmark (PR #110) | OpenAMP CI |
