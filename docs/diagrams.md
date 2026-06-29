# OpenAMP Foundry — Pipeline Diagrams & Expert Reasoning

> **Purpose.** A single source of truth for *how* the dry-lab pipeline turns random
> sequence space into a small set of synthesis-ready, novel, selective antimicrobial
> peptide (AMP) candidates — and *why* each stage exists. Written to be defensible to a
> peptide chemist, a computational biologist, and an IP attorney simultaneously.
>
> **Honest framing.** Every score here is a *transparent heuristic*, not a validated
> biological predictor. The pipeline's job is to **narrow 10⁷ → ~10¹ with an auditable
> trail**, not to claim discovery. Discovery requires wet-lab confirmation. Nothing in
> this document asserts a candidate *works* — only that it survives every cheap filter
> a domain expert would apply before spending money on synthesis.

---

## 1. The Core Problem (why this is hard)

```
                  Sequence space for a 12–24-mer peptide
                  ≈ 20^18  ≈ 10^23 possible sequences
                              │
            ┌─────────────────┴──────────────────┐
            │  Almost all are useless:            │
            │   • not antimicrobial               │
            │   • or antimicrobial BUT hemolytic  │  ← the real trap
            │   • or active BUT already known     │  ← novelty/IP trap
            │   • or novel BUT unsynthesisable    │
            └─────────────────┬──────────────────┘
                              ▼
         We must find the vanishingly rare sequences that are
     ACTIVE  ∧  SELECTIVE  ∧  NOVEL  ∧  SYNTHESISABLE  ∧  IP-CLEAR
```

**Expert insight that shapes everything below:** the binding constraint is **not**
"is it an AMP" (predictors say yes easily) — it is **"will it spare host cells"**
(selectivity / low hemolysis). So the pipeline is deliberately weighted toward
selectivity and safety, and rejects the melittin-like high-hydrophobic-moment helices
that every naive generator over-produces.

---

## 2. End-to-End Pipeline (bird's-eye)

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 0   KNOWLEDGE BASE          51,503 unique known AMPs (6 databases)        │
│           build_db() + build_kmer_index()                                       │
└───────────────┬──────────────────────────────────────────────────────────────┘
                │  (loaded once; reused by novelty scan AND motif prior-art)
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 1   GENERATION              biased random sampling from weighted alphabet │
│           generate_expert_1000.py   ~10^5–10^6 raw sequences                     │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 2   EXPERT GATES (cheap, main process)        cost-ordered, reject early │
│   2a biophysics   charge, μH band, hydrophobic+aromatic caps, Trp≤2             │
│   2b selectivity  selectivity_proxy ≥ 0.60   (charge/GRAVY therapeutic window)  │
│   2c synthesis QC  no DKP / aspartimide / Trp-photolability; difficulty ≠ HIGH  │
│   2d motif novelty k-mer prior-art: no ≥3 consecutive known 5-mers              │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼  (~15–18% survive)
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 3   NOVELTY SCAN (expensive, 9 parallel workers)                          │
│           BLOSUM62 local identity vs all 51,503 AMPs;  keep < 40%               │
│           patent DB flagged → any DRAMP-patent proximity rejected (CLEAR only)  │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 4   EXPERT COMPOSITE + DIVERSITY                                          │
│           expert_score(): activity∩selectivity∩safety∩synth∩motif∩hinge        │
│           diversity bucketing: ≤28 per (charge × length × hydrophobicity) bin   │
│           → 1000 ranked, scaffold-diverse candidates                            │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 5   EXTERNAL VALIDATION (independent models we do NOT control)            │
│   local : Macrel (AMP + hemolysis)                                              │
│   web   : AMPScanner v2 · CAMPR4 · HemoFinder · AntiCP2                          │
│           strict shortlist = AMP+ ∧ NonHemo ∧ low-AntiCP ∧ predictor-consensus  │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 6   IP CLEARANCE (pre-filing)                                             │
│   BLASTp vs NCBI nr + pataa  ·  Lens.org / WIPO PatentScope exact-string        │
└───────────────┬──────────────────────────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ STAGE 7   WET-LAB PANEL          8–10 peptides → MIC + hemolysis + cytotoxicity │
│           the ONLY stage that can confirm discovery                             │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Why Cost-Ordering (the single most important engineering choice)

A BLOSUM62 scan against 51,503 sequences costs ~10–100 ms per candidate; computing
biophysical features + QC + k-mer lookup costs <1 ms. If we scanned every raw sequence
we would waste >90% of compute on candidates that fail a trivial charge or motif check.

```
  Raw sequence
      │
      ├─►  cheap gates  (µs–ms)  ──►  reject ~82–85%   ◄── do this FIRST
      │
      └─►  survivors only ──►  BLOSUM62 scan (ms, 9 workers)  ◄── expensive, LAST
```

This is why Stage 2 lives in the main process and only Stage 3 is parallelised.
**Reject the cheap-to-reject things first.**

---

## 4. The Expert Composite (Stage 4 detail)

No single number captures a good candidate. A 30-year peptide chemist holds several
axes at once; `scoring/expert.py` makes that explicit and auditable.

```
expert_composite =
        0.22 · activity_consensus     physchem ∩ Boman, minus disagreement penalty
      + 0.22 · selectivity            charge/GRAVY therapeutic-window proxy
      + 0.18 · safety                 hemolysis-risk proxy (μH, hydrophobicity, charge)
      + 0.13 · synthesis              SPPS feasibility (length, repeats, aggregation)
      + 0.05 · serum_stability        proteolytic longevity (informational)
      + 0.08 · hinge_selectivity      central Gly/Pro hinge motif  ◄── NEW
      + 0.12 · motif_novelty          k-mer prior-art beyond global identity ◄── NEW
        ─────
         1.00   (weights fixed BEFORE ranking; deliberately tilted to selectivity+safety)
```

### 4a. Why two NEW signals that no prior module had

**Helix-hinge selectivity** — *"the expert sees a break at position N."*

```
   Rigid amphipathic helix              Hinged helix (Gly/Pro in centre)
   ───────────────────────              ────────────────────────────────
   ███████████████████████              ██████████▼██████████
   continuous hydrophobic face          two short faces, flexible kink
        │                                     │
        ▼                                     ▼
   efficient at lysing BOTH              still attracted to anionic
   bacterial AND zwitterionic            bacterial membranes, but a poor
   mammalian membranes → HEMOLYTIC       continuous pore-former in
                                         zwitterionic mammalian membranes
                                         → BETTER THERAPEUTIC WINDOW
```

A single helix-breaker (G/P) in the central third (cecropin-A's Gly23-Pro24 class)
scores 1.0; a continuous run of breakers shreds the helix and scores 0.1. A flat
hydrophobic-moment number is **blind to position** — this feature is not.
*Refs: Tossi et al. 2000 Biopolymers; Shai 2002 Biopolymers; Saberwal & Nagaraj 1994 BBA.*

**Motif-level prior art** — *"the expert recognises the local sequence."*

```
   Global BLOSUM identity: 38%   →  "looks novel"  ✔ passes the identity gate
                    BUT
   contains  ...R-W-W-K-G-G-W-Q...  a 7-mer copied verbatim from a known AMP
                    →  k-mer prior-art flags it as locally derivative  ✘
```

`build_kmer_index()` stores every 5-mer in the 51,503-sequence corpus; a candidate
sharing ≥3 *consecutive* known 5-mers is rejected. This is strictly sharper than a
global-identity threshold for catching "this looks like X."

### 4b. The whole per-axis toolkit feeding the composite

```
   compute_features(seq) ─┬─► activity_likeness_score   (Zasloff, Hancock&Sahl)
                          ├─► boman_activity_score       (Boman 2003)  ← independent
                          ├─► model_disagreement         (uncertainty penalty)
                          ├─► safety_score               (Dathe&Wieprecht 1999; μH)
                          ├─► synthesis_feasibility      (SPPS rules)
                          ├─► serum_stability_score      (Hilpert 2006 protease sites)
                          ├─► selectivity_proxy          (Shai 2002 charge/GRAVY)
                          ├─► helix_hinge_analysis        ◄── NEW
                          └─► kmer_prior_art              ◄── NEW
   check_sequence(seq) ──► pre-synthesis QC: aspartimide(DG/DS), deamidation(NG/NS),
                           DKP(X-Pro), pyroglutamate(Q1), Trp-photolability(≥3W),
                           aggregation runs, C-amidation & N-acetylation guidance,
                           proline-rich intracellular-assay flag
```

---

## 5. Novelty & IP — defence in depth

```
  Layer 1  k-mer prior-art      local motif lifted from known AMP?      (generation)
  Layer 2  BLOSUM62 vs 51,503   global similarity < 40%?  patent-flagged? (generation)
  Layer 3  BLASTp NCBI nr+pataa  global vs ALL proteins + patent proteins (pre-filing)
  Layer 4  Lens.org / WIPO       exact-string in patent claims/listings   (pre-filing)
```

Thresholds (identity = matches / query length, BLOSUM62 local):

```
  ≥99%  EXACT_MATCH_OR_FRAGMENT     ┐
  ≥80%  KNOWN_VARIANT               │ excluded from novelty claims
  ≥60%  CLOSE_RELATIVE              ┘
  ≥40%  RELATED_NOVEL               ← allowed; note similarity in docs
  <40%  HIGH_CONFIDENCE_NOVEL       ← target
  patent DB hit at ANY identity ≥60% → POSSIBLE_PATENT_RISK → SYNTHESIS_HOLD
```

---

## 6. What each stage can and cannot claim

| Stage | Can claim | Cannot claim |
|-------|-----------|--------------|
| 1–2 Generation+gates | "biophysically AMP-plausible, synthesisable" | any activity |
| 3 Novelty | "<40% identity to 51,503 known AMPs; no patent proximity" | global IP clearance |
| 4 Expert composite | "ranks well on transparent multi-axis heuristic" | biological potency |
| 5 External predictors | "independent models concur it's AMP-like and low-hemolysis" | in-vitro activity |
| 6 IP clearance | "no exact prior art found in searched databases" | freedom-to-operate (needs counsel) |
| 7 Wet lab | **"experimentally validated"** | (this is the real result) |

> **The hard ceiling.** Every dry-lab improvement raises *confidence and efficiency*,
> not *proof*. The unlock is Stage 7. The pipeline exists to make Stage 7 cheap, fast,
> and aimed at the most defensible candidates — and to publish the negative results,
> benchmarks, and leakage checks regardless of outcome.

---

## 7. Data & Code Map

| Stage | Code | Output |
|-------|------|--------|
| 0 | `build_db()`, `build_kmer_index()` in `scripts/generate_expert_1000.py`; DBs in `data/novelty_db/` | in-memory corpus + 5-mer index |
| 1–4 | `scripts/generate_expert_1000.py` · `src/openamp_foundry/scoring/expert.py` | `outputs/expert_1000_candidates.{csv,fasta}` |
| 2c | `src/openamp_foundry/qc/presynth_check.py` | synthesis liability flags |
| 3 | `scripts/run_expanded_novelty_audit.py` | independent novelty re-audit |
| 5 | `scripts/screen_1000_candidates.py` (Macrel) + web predictors | `outputs/screening_1000/` |
| 6 | `scripts/run_patent_blastp.py` | `outputs/patent_blastp_*/ip_clearance_report.md` |
| 7 | `docs/WET_LAB_HANDOFF.md` | synthesis order + assay plan |

*See `docs/NOVELTY_AUDIT_GUIDE.md` for the canonical novelty methodology and database provenance.*
