# PRKCB project v2 — build plan

**Goal:** a pipeline that demonstrates professional bioinformatics practice, plus
one genuine methodological finding.

**Time:** Thu 19 Nov (part), Fri 20, Wed 25, Thu 26, Fri 27.

**Rule:** nothing gets written until you've approved it. Each day ends with
something that runs and something you can explain.

---

# Part 1 — What "professional" actually means here

The v1 project failed a specific test. Not correctness — the answers were right —
but **method**. Three things in it were done by hand that the field does with
standard tools:

| v1 did this | The field does this | Why it matters |
|---|---|---|
| Homemade regex to find domains | HMMER / InterProScan | Regex has no statistics. No E-value, no null model |
| Thermodynamic values typed from memory | ViennaRNA | Hand-typed numbers can be silently wrong |
| Substring search for off-targets | BLAST | Substring search finds only perfect matches |

A reviewer who sees regex-based domain annotation stops reading. Not because the
answer is wrong, but because it signals someone who doesn't know the standard
tools exist. That is the single highest-value thing to fix.

**Success criteria for v2.** At the end you should be able to say, truthfully:

1. Every annotation comes from a validated, citable source or tool.
2. Every thermodynamic number comes from a published parameter set via a library.
3. Off-target screening used BLAST, not string matching.
4. The pipeline found something the standard tools don't check for.
5. You can explain any line of it.

---

# Part 2 — The tools, and why each one

This section exists because using a tool you can't justify is as bad as not
using it.

## 2.1 HMMER / InterProScan — domain annotation

**The problem it solves.** Given a protein sequence, which known domains does it
contain and where?

**Why not regex.** A regex is a rigid pattern: it matches or it doesn't. Real
domains tolerate substitution — a C1 zinc finger in one protein differs from
another's at most positions while keeping its fold. More importantly, a regex
gives you **no measure of surprise**. If I search a 671-residue protein for a
loose pattern, a hit may mean nothing; a regex can't tell you that.

**What a profile HMM does instead.** A **Hidden Markov Model** profile is built
from an alignment of thousands of homologous domains. For *every position* in
the domain it stores a probability distribution over all 20 amino acids, plus
probabilities for insertions and deletions. Scoring a sequence against it asks:
*how much more likely is this sequence under the domain model than under a
random-sequence model?* That ratio becomes a bit score and an **E-value** — the
number of hits this good you'd expect by chance in a database this size.

An E-value of 1e-40 means "this is not a coincidence." A regex can never say
that.

**Pfam** is the curated database of these profiles (~20,000 domain families).
**HMMER** is the program that scores sequences against them (`hmmscan`).
**InterProScan** runs HMMER against Pfam *plus* about thirteen other member
databases — PROSITE, SMART, CDD, PRINTS, PANTHER — and integrates the results
into InterPro entries. It's the standard "annotate this protein" answer.

**What we'll do.** You run InterProScan once on the EBI website and download the
TSV. We parse it. That gives curated, citable domain boundaries with E-values.

**If that fails** (site down, output unusable): fall back to UniProt P05771's
curated feature table, which is expert-reviewed annotation. Weaker — it's a
lookup rather than a computation — but still legitimate and citable.

**Honest limitation.** For a well-studied human kinase this is a *solved*
problem. We are not discovering domains. We are using the right tool to get the
right answer, and being able to say so.

**`pyhmmer`** (already installed, works) lets us run HMMER from Python if we
later want to score sequences against specific Pfam profiles directly. It needs
the Pfam HMM files, which you'd have to download — only worth it if InterProScan
isn't enough.

## 2.2 ViennaRNA — RNA structure and thermodynamics

Already installed and tested (v2.8.0). This one replaces two weaknesses at once.

**Problem 1 — the thermodynamic parameters.** v1 contained a table of RNA
nearest-neighbour free energies typed in by hand. ViennaRNA ships the full
**Turner 2004** parameter set: not just the 16 nearest-neighbour stacks but
hairpin loops, bulges, internal loops, dangling ends, terminal mismatches —
hundreds of experimentally measured values. Using the library means the numbers
are right and you can cite them.

**Problem 2 — accessibility, which v1 ignored completely.**

An mRNA is not a straight line. It folds into hairpins and loops. A target site
buried inside a stable double-stranded stem is **physically inaccessible** —
RISC can't reach it — no matter how good its composition score. v1 had no
concept of this at all.

ViennaRNA computes:

- **MFE (minimum free energy) structure** — the most stable fold, in dot-bracket
  notation where `.` is unpaired and `(` `)` are paired bases.
- **Base-pair probabilities** via the partition function — better than MFE
  alone, because it accounts for the ensemble of structures rather than one.
- **Accessibility / opening energy** — the energy cost of forcing a region to be
  unpaired, which is what a reagent must pay to bind.
- **Duplex energy** (`RNA.duplexfold`) — the energy of guide binding to target,
  which replaces the hand-rolled calculation entirely.

**How we'll use it.** For each candidate, fold a window of transcript around the
target site and compute how accessible that site is. Candidates in stable stems
get penalised. This adds a dimension v1 simply didn't have.

**If it's too slow** (folding thousands of windows): fold a local window
(~80 nt) rather than the whole transcript, which is standard practice anyway
since local structure dominates, and cache results.

**Honest limitation.** Computational folding predicts *in vitro* structure. Real
mRNA in a cell is coated in proteins and being translated, which changes
everything. Accessibility scores are a useful filter, not truth.

## 2.3 BLAST — off-target screening

**The problem.** v1 asked "does this exact 19-mer appear elsewhere?" using
Python's `in` operator. That finds only **perfect matches**. Real off-target
binding tolerates mismatches — especially outside the seed. So v1's off-target
screening was close to meaningless, which the README admits.

**What BLAST does.** Basic Local Alignment Search Tool finds *similar*
subsequences, not identical ones. It seeds on short exact words, extends them
into local alignments, scores with a substitution model, and reports E-values.
That's what you need to ask "what else in the transcriptome might this bind,
allowing for mismatches?"

**For short queries like a 19-mer**, standard BLAST settings miss things — the
default word size is too large. The right settings are `blastn-short`
(`-task blastn-short`), a small word size (7), and a relaxed E-value threshold.
Knowing that is itself part of using the tool properly.

**The practical problem.** A proper screen needs the human transcriptome
database (~200 MB+), and this session cannot download it.

**Two routes, and I'd like your view:**

- **(a) Web BLAST.** You run the top ~20 candidates against human RefSeq RNA on
  the NCBI website and save the hit tables. Real screening against the real
  database. Costs you maybe an hour of clicking. Limited to the shortlist.
- **(b) Local BLAST against a curated set.** You download ~15 transcripts (PKC
  family and close relatives); I build a small BLAST database locally and screen
  *all* candidates against it. Comprehensive within a narrow universe.

Best answer is **both**: (b) for systematic screening of everything, (a) for the
final shortlist. That's also exactly what a real project would do — cheap broad
filter, expensive precise check on the finalists.

## 2.4 Supporting tools

**Biopython** — sequence parsing, translation, and `MeltingTemp` for ASO hybrid
Tm using the published Sugimoto table. Already in use, already appropriate.

**pandas** — for the candidate tables. v1 used `csv.DictWriter` and manual
loops; pandas is what the field actually uses for tabular work, makes filtering
and grouping readable, and signals familiarity with the standard stack.

**pytest** — v1's tests were hand-rolled `check()` calls. pytest is the standard
framework. Same tests, recognisable structure, and `pytest -v` output is
something a reviewer immediately understands.

**git** — proper commit history this time. Small commits with real messages, not
one "Add files via upload."

---

# Part 3 — The depth piece, specified

*Revised 19 Nov after literature validation. The first version of this section
was built on a wrong mechanism; see 3.7.*

## 3.1 What the literature establishes

Two separate mechanisms, which v1 of this plan conflated:

**Cleavage discrimination.** When an siRNA must distinguish two near-identical
targets, the position of the mismatch decides whether it succeeds. Huang et al.
(2009, *NAR* 37:7560) profiled this position by position and found the
**cleavage-site region, roughly guide positions 9–13**, is most sensitive —
position 12 especially. They report "generally low sensitivity was found in the
seed region," and position 17 is also low-sensitivity. So the seed is
comparatively *tolerant* of mismatches for this purpose.

**Seed-mediated off-target repression.** Separately, perfect seed complementarity
(guide positions 2–8) causes miRNA-like repression of unintended transcripts
(Birmingham et al. 2006, *RNA* 12:1179). Different mechanism, different question.

Both matter here, but they are not the same thing and must be scored separately.

## 3.2 The geometry specific to this gene

Because βII's 3′ UTR retains the entire βI exon, a junction-spanning siRNA has
**two** potential off-target sites in βII, not one:

1. βII's **shared exon** — matches the shared portion of the window
2. βII's **retained βI exon** in the 3′ UTR — matches the βI portion

With the junction at window position *j*, the longest contiguous match to βII is
therefore `max(j, 19−j)`.

Since `junction_balance = min(j, 19−j)`, this gives the exact identity:

> **longest contiguous match to βII = 19 − junction_balance**

Verified across all 18 candidates with no exceptions.

**This is the part that standard allele-specific design rules do not cover.**
They assume a single off-target site. Here the competing transcript contains
*both* flanks of the junction, so discrimination requires minimising the longer
flank — which means placing the junction near the middle of the window.

## 3.3 The hypothesis

> A junction-spanning siRNA discriminates only when the junction sits near the
> centre of the window, so that neither flank is long enough to support
> productive binding to the competing transcript. Standard composition scores
> are blind to this and will rank non-discriminating candidates highest.

## 3.4 Verified against v1 data, 19 Nov

| Candidate | junction_balance | Longest βII match | Composite score |
|---|---|---|---|
| beta1_2054 — v1's top pick | 3 | **16** | 75.5 |
| beta1_2052 | 5 | 14 | 66.1 |
| beta1_2048 — actually best | 9 | **10** | 63.6 |
| beta1_2047 — actually best | 9 | **10** | 46.4 |

The identity holds exactly for all 18. The ranking is inverted: v1's headline
recommendation has the *longest* match to the competing transcript of any
high-scoring candidate.

## 3.5 What v2 computes

1. Junction position *j* and `junction_balance` per candidate.
2. Longest contiguous match to βII, at **both** off-target sites separately.
3. **Position-weighted mismatch score** using Huang et al.'s sensitivity profile
   — weight mismatches at guide 9–13 heavily, seed and position 17 lightly.
   This replaces the naive "count the mismatches" approach.
4. Seed-match screening against βII, scored as a **separate** liability
   (miRNA-like repression, not failed cleavage).
5. Mismatch-tolerant off-target search by BLAST, not substring matching.
6. Rank correlation between composite score and position-weighted discrimination.

## 3.6 Possible outcomes, named in advance

- **Standard scoring is inverted** — expected, already visible. Produces a
  concrete recommendation: rank junction candidates by balance first.
- **No βI candidate is usable** — a real possibility. The two best by geometry
  score 63.6 and 46.4 on composition and still need to pass accessibility. If
  none survives, "βI cannot be selectively targeted by siRNA; use a gapmer or
  accept pan-β knockdown" is a clean, honest result.
- **Something passes everything** — then there is a defensible lead compound.

All three are reportable. None is a failed project.

## 3.7 Novelty, stated honestly

Position-dependent mismatch discrimination is **established** (Huang 2009;
Schwarz et al. 2006, *PLoS Genet* 2:e140). "Where the mismatch sits matters" is
not a new idea and must not be presented as one.

What is defensible: applying those established principles to an
alternative-last-exon gene where the competing transcript retains both flanks,
which creates a **two-sided constraint** that single-off-target design rules do
not address. `junction_balance` is the metric that falls out of that geometry.

The claim is: *correct application of known principles to an unusual splicing
geometry, where the standard tools give actively misleading answers.* Small, but
true and demonstrable.

## 3.8 Extension, only if time allows

Repeat on 3–5 other alternative-last-exon genes. First thing cut.

# Part 4 — Day by day

Each day has an objective, an acceptance test, and a contingency.

## Today (Thu 19) — setup

**Objective:** folder skeleton, data manifest, and you start downloading.

**Build:** directory structure, `requirements.txt` with pinned versions,
`data/MANIFEST.md` recording exactly what each file is and where it came from,
empty git repo with a real first commit.

**Acceptance:** `pip install -r requirements.txt` works from clean; manifest
lists every file we will need.

**If short on time:** manifest only. Everything else is quick.

## Fri 20 — the annotation layer

**Objective:** kill the regex. Replace with validated annotation.

**Needs from you:** InterProScan TSV, UniProt P05771 entry.

**Build:** a parser for the InterProScan output; a rewritten script 01 that
reports domains with their source database, accession and E-value; a comparison
table showing where v1's regex coordinates agreed and disagreed with the curated
ones.

**Why the comparison table matters:** it turns an embarrassment into a result.
"Here is what naive pattern matching gets wrong relative to profile HMMs" is a
legitimate methods observation, and it's honest about the project's history.

**Acceptance:** every domain coordinate traceable to an InterPro accession.

**If InterProScan output doesn't arrive:** fall back to the UniProt feature
table. Note the weaker provenance explicitly.

## Wed 25 — structure and thermodynamics

**Objective:** real energetics and accessibility.

**Build:** replace the hand-typed `NN_DG` with ViennaRNA throughout; add
accessibility scoring per candidate; re-score everything; compare old and new
rankings.

**Acceptance:** no hand-entered thermodynamic constant anywhere; every candidate
carries an accessibility value.

**Contingency:** if folding every candidate is too slow, restrict accessibility
to candidates that already pass composition filters. Document the restriction.

**Watch for:** rankings will change. That's the point — but we must be able to
explain *why* each change happened, not just note it.

## Thu 26 — off-target and the depth piece

**Objective:** BLAST screening, then the junction-geometry analysis.

**Needs from you:** paralog transcripts; optionally web BLAST results for the
shortlist.

**Build:** local BLAST database and screening of all candidates; the junction-geometry
analysis from Part 3; the figure that shows it.

**Acceptance:** every candidate has a BLAST-based off-target count; the
junction-geometry hypothesis is tested and the answer is stated either way.

**This is the day with the most risk.** If BLAST setup eats the morning, the
depth piece moves to Friday and figures get compressed. The depth piece does
**not** get cut — it's the reason we chose this plan.

## Fri 27 — write-up and defence

**Objective:** make it presentable and make you able to defend it.

**Build:** figures, README rewritten around what v2 actually did, a methods
section with every tool, version and citation, and clean git history.

**Then:** I ask you hard questions about your own project until you can answer
them without looking.

---

# Part 5 — Data manifest

| # | What | Where | Needed by | Purpose |
|---|---|---|---|---|
| 1 | InterProScan TSV for PKC-βI protein | ebi.ac.uk/interpro/search/sequence/ | **Fri 20** | Domain annotation |
| 2 | UniProt P05771 full entry (text) | uniprot.org | **Fri 20** | Curated features, fallback |
| 3 | PRKCA NM_002737 transcript FASTA | NCBI Nucleotide | Wed 25 | Off-target |
| 4 | PRKCG NM_002739 transcript FASTA | NCBI Nucleotide | Wed 25 | Off-target |
| 5 | Mouse Prkcb NM_008855 transcript | NCBI Nucleotide | Wed 25 | Cross-species, nucleotide-level |
| 6 | PRKCD, PRKCE, PRKCZ transcripts | NCBI Nucleotide | Thu 26 | Wider off-target set |
| 7 | Web BLAST hit tables for shortlist | blast.ncbi.nlm.nih.gov | Thu 26 | Real transcriptome screen |

Items 1 and 2 are the critical path. Everything else has slack.

Already in hand: βI and βII proteins, mouse protein, both PRKCB transcripts.

---

# Part 6 — Triage, in cut order

If time runs short, things get dropped in this order. Decided now, while
nobody's panicking:

1. Extension of the depth piece to other genes *(optional from the start)*
2. Items 6–7 of the manifest — wider off-target set
3. ASO arm — keep siRNA, note ASO as future work
4. Accessibility scoring — keep ViennaRNA for thermodynamics only
5. The v1-vs-v2 annotation comparison table

**Never cut:** validated domain annotation, ViennaRNA thermodynamics, the
junction-geometry analysis, honest limitations. Those four are the project.

---

# Part 7 — What you'll be able to say at the end

> "I built a pipeline that designs isoform-selective RNA reagents against PRKCB.
> Domains are annotated with InterProScan rather than hand-written patterns.
> Thermodynamics come from ViennaRNA's Turner parameters, and I score target-site
> accessibility by local folding. Off-target screening is BLAST-based with
> short-query settings, not substring matching.
>
> The point of it is about junction geometry. PRKCB's βII transcript retains the
> entire βI terminal exon in its 3' UTR, so an siRNA spanning the splice junction
> has two possible off-target sites in βII, one on each side of the junction. Its
> longest contiguous match to the competing transcript is 19 minus how centred the
> junction is — an identity that holds for every candidate. Discrimination
> therefore needs the junction near the middle of the window, and standard
> composition scoring is blind to that: it ranked the candidate with the *longest*
> match to the competing transcript highest.
>
> Position-dependent mismatch discrimination itself is established — Huang 2009,
> Schwarz 2006. What's mine is applying it to a splicing geometry where the
> competing transcript carries both flanks, which single-off-target design rules
> don't cover. It's computational and untested, and the splicing mechanism was
> described in 2004."

Every clause of that is something you can back with a citation or a line of
output.
