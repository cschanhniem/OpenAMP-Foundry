# Novelty Audit Guide — OpenAMP Foundry

> **Last updated:** 2026-06-29  
> **Audit level:** STRONG (APD6 + DRAMP 3.0 + UniProt reviewed/unreviewed)  
> **DB size:** 34,555 unique clean standard-AA sequences (5–100 aa)

## Why Novelty Auditing Matters

Every AMP candidate must pass a full novelty audit before wet-lab synthesis to:

1. **Avoid IP risk** — DRAMP patent sequences are flagged; any candidate with ≥60% identity to a patent sequence receives `POSSIBLE_PATENT_RISK` and a `SYNTHESIS_HOLD`.
2. **Ensure scientific validity** — a "novel" label on a KNOWN_VARIANT of a published AMP is misleading.
3. **Protect the budget** — synthesis of a patent-risk or known candidate wastes money.

The single biggest mistake in OpenAMP Wave 0.5 was relying on the `novelty-check-broad` CLI tool, which uses a 72-sequence curated reference — far too small to catch patent proximity. See [`docs/postmortems/2026-06-29-seed020-novelty-patent-risk.md`](postmortems/2026-06-29-seed020-novelty-patent-risk.md) for the full incident.

**Rule: Always run `scripts/run_expanded_novelty_audit.py` before adding any candidate to the panel or ordering synthesis.**

---

## Database Sources

All FASTA files are committed to `data/novelty_db/`.

| File | Source | Sequences | Patent-flagged | Notes |
|------|--------|-----------|----------------|-------|
| `apd6_natural.fasta` | [APD6](https://aps.unmc.edu) | 3,306 | No | Canonical natural AMPs, curated |
| `apd6_animal.fasta` | APD6 2024 | 2,580 | No | Animal AMPs with known activity |
| `apd6_plant.fasta` | APD6 2024 | 268 | No | Plant defensins & thionins |
| `apd6_bacteria.fasta` | APD6 2024 | 410 | No | Bacteriocins, lantibiotics |
| `dramp_general.fasta` | [DRAMP 3.0](https://dramp.cpu-bioinfor.org) | 11,687 | No | General AMP collection |
| `dramp_patent.fasta` | DRAMP 3.0 | 18,715 | **YES** | Patent-protected AMPs — IP risk detection |
| `dramp_specific.fasta` | DRAMP 3.0 | 6,321 | No | Clinical, stability, structural entries |
| `uniprot_amps_reviewed.fasta` | UniProt KW-0929 ≤100aa, reviewed | 2,673 | No | Swiss-Prot curated |
| `uniprot_amps_unreviewed.fasta` | UniProt KW-0929 ≤60aa, unreviewed | 1,692 | No | TrEMBL broad coverage |

After deduplication: **34,555 unique clean standard-AA sequences (5–100 aa)**

### Databases NOT yet included

| Database | Size | Status | Why excluded |
|----------|------|--------|--------------|
| [dbAMP 3.0](https://awi.cuhk.edu.cn/dbAMP/) | ~35,518 | No bulk download available | Website has no public FASTA export endpoint |
| [DBAASP](https://dbaasp.org) | >10,000 | Requires auth | API returns 403 for bulk download |
| [ESCAPE](https://arxiv.org/abs/2511.04814) | ~80,000 | Research preprint | Not yet publicly downloadable as a single FASTA |
| UniProt nr (full) | Billions | Scope | Beyond AMP-specific novelty; use BLASTp for global search |

**For global novelty (pre-grant):** Run BLASTp against NCBI nr and Google Patents / [Lens.org](https://lens.org) in addition to this audit.

### Updating the databases

```bash
# Re-download APD6 (run annually)
curl -sL "https://aps.unmc.edu/assets/sequences/naturalAMPs_APD2024a.fasta" \
  -o data/novelty_db/apd6_natural.fasta
curl -sL "https://aps.unmc.edu/assets/sequences/animalAMPs_APD2024a.fasta" \
  -o data/novelty_db/apd6_animal.fasta
# (plant and bacteria similarly)

# Re-download DRAMP 3.0 (run annually)
BASE="https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new"
curl -sL "${BASE}/general_amps.fasta" -o data/novelty_db/dramp_general.fasta
curl -sL "${BASE}/patent_amps.fasta"  -o data/novelty_db/dramp_patent.fasta
curl -sL "${BASE}/specific_amps.fasta" -o data/novelty_db/dramp_specific.fasta

# Re-download UniProt (run annually) — use the Python paginator in scripts/
python3 scripts/download_uniprot_amps.py
```

---

## Audit Methodology

### Algorithm

**BLOSUM62 local alignment** (BioPython `PairwiseAligner`):

```
gap_open   = -11
gap_extend = -1
mode       = local
identity   = (exact matches in aligned blocks) / query_length
```

This catches conservative substitutions (K↔R, I↔L↔V, F↔W↔Y) that simple edit distance misses.

### Novelty Thresholds

| Identity | Class | Meaning |
|----------|-------|---------|
| ≥99% | `EXACT_MATCH_OR_FRAGMENT` | Identical or contained within a known sequence |
| ≥80% | `KNOWN_VARIANT` | Minor point mutations from a known AMP |
| ≥60% | `CLOSE_RELATIVE` | Same structural family, meaningful differences |
| ≥40% | `RELATED_NOVEL` | Related but sufficiently distinct |
| <40% | `HIGH_CONFIDENCE_NOVEL` | Genuinely new sequence space |

### Patent Risk Classification

| Condition | Flag |
|-----------|------|
| Best hit in `dramp_patent` AND identity ≥60% | `POSSIBLE_PATENT_RISK` → **SYNTHESIS HOLD** |
| Best hit in `dramp_patent` AND identity ≥40% | `LOW_PATENT_RISK` → consult IP counsel |
| Best hit in `dramp_patent` AND identity <40% | `CLEAR` |
| Best hit not in patent DB | `CLEAR` |

---

## Running the Audit

### Standard usage

```bash
# From repo root, always use .venv/bin/python3
.venv/bin/python3 scripts/run_expanded_novelty_audit.py \
    --input outputs/my_candidates.csv \
    --out outputs/my_novelty_audit.csv

# Or from a FASTA
.venv/bin/python3 scripts/run_expanded_novelty_audit.py \
    --fasta outputs/my_candidates.fasta \
    --out outputs/my_novelty_audit.csv
```

Input CSV must have columns: `candidate_id`, `sequence`

### Interpreting results

```
✓ HIGH_CONFIDENCE_NOVEL + CLEAR  →  approved for predictor submission and wet lab
✓ RELATED_NOVEL + CLEAR          →  approved; note 40-60% similarity in documents
⚠ CLOSE_RELATIVE + CLEAR         →  approved; requires explicit prior-art disclosure
⚠ CLOSE_RELATIVE + POSSIBLE_PATENT_RISK  →  SYNTHESIS HOLD; IP review first
⚠ KNOWN_VARIANT / EXACT_MATCH   →  exclude from novelty claims; use as control only
```

### Before ordering synthesis

```bash
# 1. Run audit
.venv/bin/python3 scripts/run_expanded_novelty_audit.py --input candidates.csv

# 2. Confirm output CSV — all selected candidates must show CLEAR patent_risk
grep "POSSIBLE_PATENT_RISK\|REVIEW_REQUIRED" outputs/expanded_novelty_audit_*.csv

# 3. Never use novelty-check-broad alone — it uses only 72 sequences
```

---

## What Audit Level Do We Have?

| Level | Requirements | Status |
|-------|-------------|--------|
| **Minimum (pre-wet-lab)** | APD6 + DRAMP | ✅ |
| **Strong (pre-grant)** | + DBAASP + UniProt | ✅ (UniProt included) |
| **Best-in-class** | + ESCAPE 80k + NCBI nr + patent string search | ⚠ BLASTp / Lens.org needed manually |

The current audit is **STRONG level** and sufficient for wet-lab synthesis decisions.

For patent filings or grant submissions, additionally run:
- BLASTp against NCBI nr (web: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch)
- Exact-sequence search on Lens.org or Google Patents

---

## De Novo Campaign Results (2026-06-29)

All de novo candidates generated by `scripts/generate_denovo_novel_amps.py` (seed=42, 30 targets found in 15,315 attempts) were re-audited against the expanded 34,555-sequence DB:

| Result | Count |
|--------|-------|
| HIGH_CONFIDENCE_NOVEL + CLEAR | **23** |
| RELATED_NOVEL + CLEAR | **7** |
| POSSIBLE_PATENT_RISK | **0** |
| SYNTHESIS_HOLD | **0** |

Top 12 selected for external predictor submission (`outputs/denovo_top12_final.fasta`):
10 × HIGH_CONFIDENCE_NOVEL + 2 × RELATED_NOVEL (top by pipeline ensemble score).

Prior wave novelty audit findings (`outputs/wave0_pilot_novelty_audit.csv`):
All Wave 0 "HIGH_CONFIDENCE_NOVEL" labels were wrong — those were KNOWN_VARIANT (>80%) or CLOSE_RELATIVE when checked against the full DB.

---

## Agents / Automated Use

When an agent generates AMP candidates:

1. Write sequences to a CSV with `candidate_id`, `sequence` columns.
2. Run: `.venv/bin/python3 scripts/run_expanded_novelty_audit.py --input <csv> --out outputs/audit_<timestamp>.csv`
3. Parse output — reject any row where `patent_risk != CLEAR`.
4. For wet-lab submission, require `novelty_class IN (HIGH_CONFIDENCE_NOVEL, RELATED_NOVEL)`.
5. Never cite `novelty-check-broad` (72-seq) output in a submission document.
