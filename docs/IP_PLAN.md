# IP and Disclosure Plan

> **Purpose:** Strategic guidance for intellectual property decisions before public
> sequence disclosure. Not a legal opinion — consult a patent attorney before filing.
>
> **Last updated:** 2026-06-29 (Sprint 6)

---

## Patent Strategy

### Claimable Candidates

| Seed | Classification | Claim strategy |
|------|---------------|----------------|
| SEED-006 (4 candidates) | HIGH_CONFIDENCE_NOVEL | Primary patent claims — novel wasp-venom scaffold |
| SEED-007 (1 candidate) | HIGH_CONFIDENCE_NOVEL | Primary patent claims — novel bumblebee-venom scaffold |
| SEED-008 (3 candidates) | HIGH_CONFIDENCE_NOVEL | Primary patent claims — novel Trp-rich plant scaffold |
| SEED-009 (2 candidates) | HIGH_CONFIDENCE_NOVEL | Primary patent claims — novel proline-rich scaffold |
| SEED-007 (3 candidates) | NOVEL | Moderate claims — verify mechanism distinct from bombolitin literature |
| SEED-008 (1 candidate) | NOVEL | Moderate claims — verify mechanism distinct from indolicidin |
| SEED-005 (1 candidate) | CLOSE_RELATIVE | Not individually patentable — SAR control |
| SEED-001 (1 candidate) | KNOWN_VARIANT | Not patentable — positive control |
| SEED-003 (2 candidates) | KNOWN_VARIANT | Not patentable — tachyplesin-like SAR |

### Provisional Patent

Before any public sequence disclosure, file provisional patent application covering:
1. SEED-006 scaffold and all pilot variants
2. SEED-007 scaffold and all pilot variants
3. SEED-008 scaffold and all pilot variants
4. SEED-009 scaffold and all pilot variants
5. Composition claims for any MIC ≤ 8 µg/mL candidate
6. Methods of use (antimicrobial treatment, wound healing)

### Sequence Deposit

Deposit sequences in patent filing as SEQ ID NOs. Use the WIPO ST.26 standard
(valid from July 2022) for sequence listing.

## Disclosure Plan

| Tier | Candidates | Disclosure | Timing |
|:----:|------------|------------|--------|
| 1 | SEED-006, 007, 008, 009 (14 candidates) | **Redacted** from public repo | After provisional filing |
| 2 | SEED-001, 003, 005 (6 candidates) | Public (known scaffolds) | Anytime |
| 3 | Pipeline code | Public (MIT license) | Already public |

## Freedom to Operate

**Required searches before publication:**
- [ ] APD3 BLASTp — no hit > 80% identity
- [ ] DRAMP patent section — no overlapping claims
- [ ] Google Patents / Lens.org — scaffold-level search
- [ ] Competitor sequence database (AMP-Designer, AMPGAN v3, etc.)

## Publication Strategy

1. **First publication:** Pipeline methodology + benchmark results only
   (no candidate sequences disclosed)
2. **Second publication:** Candidate nominations + novelty audit
   (after provisional patent filing)
3. **Third publication:** Wet-lab results
   (after independent replication)

## Material Transfer

Before sharing candidate peptides with external labs, execute an MTA covering:
- Sequence confidentiality
- Publication rights
- IP ownership (background vs foreground)
- Data sharing obligations

## Contributor IP

All contributors to this repository have assigned code contributions under
the MIT license. No separate contributor agreement is in place for computational
contributions. If wet-lab collaborators join, a formal collaboration agreement
should define IP ownership before any experimental work begins.
