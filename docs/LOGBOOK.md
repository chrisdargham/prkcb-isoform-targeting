# Logbook

A dated record of what was done, why, and what came out. One entry per working
session, newest at the bottom.

**What goes in an entry:** what you set out to do, what you actually did, the
key numbers, anything that surprised you, and what's next. **What doesn't:**
anything you can read off the code. This is the reasoning, not the
implementation.

**Why keep it.** Three reasons, in increasing order of how much they bite:

1. In two weeks you will not remember why a threshold is 0.01 rather than 0.05.
2. When someone asks "why did you do it that way?", this is the answer.
3. Results without a record of how they were produced are not results.

The git history records *what changed*. This records *why*, and what the numbers
were on the day.

---

## 2026-11-19 (Thu) — project restarted as v2

**Why there is a v2.** v1 reached correct conclusions by amateur means. The
specific trigger was asking how domains had been located: the answer was
hand-written regular expressions, which is circular — the patterns were derived
from the sequence they then searched. That produces no E-value, no null model,
and reports the extent of the pattern rather than the extent of the domain.

**Decided.** Rebuild with standard tools, keeping v1's analysis architecture
(which is sound) and replacing its tooling (which is not). Scope: demonstrate
professional practice, plus one methodological finding. Deadline end of
November, working Wed/Thu/Fri.

**Literature validation, before committing to the plan.** Asked whether the
planned "depth piece" was supportable. It was not, as originally framed:

- I had claimed discrimination requires the splice junction to fall inside the
  siRNA seed (guide positions 2–8).
- Huang et al. 2009 (*NAR* 37:7560) profiled mismatch discrimination by
  position and found the opposite: *"generally low sensitivity was found in the
  seed region"*. The cleavage-site region, guide positions ~9–13, dominates.
- Two mechanisms had been conflated: seed complementarity drives miRNA-like
  off-target repression (Birmingham et al. 2006), while cleavage discrimination
  is governed by the cleavage site.

**Corrected mechanism, which is better.** Because βII's 3′ UTR retains the whole
βI exon, a junction-spanning siRNA has **two** potential off-target sites in
βII — one on each side of the junction. Its longest contiguous match to βII is
`max(j, 19−j)`, where *j* is the junction position in the window. Since
`junction_balance = min(j, 19−j)`:

> **longest match to βII = 19 − junction_balance**

Verified across all 18 βI junction candidates from v1 with no exceptions.
Discrimination therefore needs the junction near the middle of the window.

`junction_balance` was the right metric in v1; the seed-based reasoning for it
was wrong. The real reason is the two-flank geometry.

**Consequence.** v1's headline βI recommendation (`beta1_2054`, composite 75.5)
has a 16-nt contiguous match to βII — near-worst available. The best by geometry
are `beta1_2047` and `beta1_2048` (balance 9, longest match 10), which score
46.4 and 63.6. Standard composition scoring is inverted relative to real
discrimination here.

**Novelty, honestly.** Position-dependent mismatch discrimination is established
(Huang 2009; Schwarz et al. 2006). What is defensible is applying it to an
alternative-last-exon gene where the competing transcript carries *both* flanks,
a two-sided constraint single-off-target rules don't address.

**Built.** Project skeleton, `requirements.txt` pinned to installed versions,
`data/MANIFEST.md` with provenance for every input.

**Next.** Annotation layer once InterProScan output arrives.

---

## 2026-11-19 (Thu, later) — annotation layer

**Inputs received.** InterProScan run on PKC-βI (671 aa): 167 hits across 13
member databases. UniProt P05771 (reviewed): 115 canonical features. Both
checksummed into `data/MANIFEST.md`.

**Built.** `src/01_annotate.py` — parses both sources, groups predictions by
InterPro entry, and reports the *spread* across databases rather than picking a
winner. Deliberate: agreement between independent methods is the evidence, and
disagreement should be visible rather than averaged away.

**Consensus boundaries.** Convergence is tight where it matters:

| Domain | Boundaries | Databases agreeing | Best E-value |
|---|---|---|---|
| C1A (phorbol/DAG 1) | 36–86 | PROSITE, Pfam, SMART | 9.4e-16 |
| C1B (phorbol/DAG 2) | 101–151 | PROSITE, Pfam, SMART | 1.9e-15 |
| C2 | 158–275 | PROSITE, Pfam, SMART | 1.2e-25 |
| Protein kinase | 342–600 (spread 3 aa) | PROSITE, Pfam, SMART | 6.0e-97 |
| AGC-kinase C-terminal | 601–671 | PROSITE, SMART | 1.3e-26 |

**The v1 comparison — the result worth keeping.**

| Feature | v1 regex | Curated | Verdict |
|---|---|---|---|
| C1A zinc finger | 37–71 | 36–86 | off by +1 / −15 |
| C1B zinc finger | 115–135 | 101–151 | off by +14 / −16 |
| C2 domain | not found | 158–275 | **missed entirely** |
| Protein kinase domain | 348–511 | 342–600 | off by +6 / −89 |
| AGC-kinase C-terminal | not found | 601–671 | **missed entirely** |
| Active site (D466) | 466 | 466 | exact |
| PDK1 site (T500) | 500 | 500 | exact |

**2 of 7 exact, both single residues. 330 residues of real domain uncovered —
49% of the protein.**

Reading: regex could find a short conserved motif but never found a *domain*. It
found where the pattern stopped. This is why the field uses profile HMMs, and
the table is now a methods observation rather than an embarrassment.

**Also built.** `src/provenance.py` — every output now carries a block
recording timestamp, git commit, whether the tree was dirty, SHA-256 of every
input, and package versions. `out/` changed from gitignored to tracked: these
files are the results and a reader should see them without re-running anything.

**Next.** ViennaRNA — replace v1's hand-typed thermodynamic table, add
target-site accessibility. Pulled forward from Wed 25 since we are ahead.

**Still needed from me:** PRKCA `NM_002737`, PRKCG `NM_002739`, mouse Prkcb
`NM_008855` transcripts.
