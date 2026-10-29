# Motivation and goal

**Status:** Standing framing doc — read alongside [phase0-spec.md](phase0-spec.md), which carries the
methodology. This doc says what counts as *finishing*.

---

## 1. The question

ichorCNA infers tumor fraction as a byproduct of fitting a copy-number mixture model, so it is blind
to near-diploid tumors; DELFI reads fragmentomic signal but only ever emits a binary cancer score.
The question this project asks sits in the gap between them:

> Can DELFI-style fragment-ratio features, trained directly against a continuous tumor-fraction
> label, recover tumor fraction where ichorCNA's CNA-based estimate is unidentifiable?

That is a falsifiable question, and the project is structured as an experiment answering it — not as
a demo of a model that works. The distinction matters for how the work is scoped, and for what
"done" means.

## 2. The goal is a conclusion, not a win

**The project succeeds if it produces a defensible answer. It does not need that answer to be yes.**

"No — fragmentomic features alone don't recover tumor fraction below ~X% in near-diploid samples,
and here is the evidence and the point at which the signal dies" is a complete result. It is a more
useful artifact than a strong R² on a model nobody can interrogate, because it demonstrates the
three things the work is actually meant to demonstrate: posing a falsifiable question, building the
machinery to answer it, and not flinching from the answer.

The in-silico admixture design is what makes either outcome defensible. Because labels are true by
construction (see [admixture-design.md](admixture-design.md)), a null result is attributable to the
*features* rather than to label noise. That is precisely what makes a negative result publishable
here rather than merely inconclusive.

## 3. What would actually sink the project

These are the real failure modes — none of them is "the model underperformed":

- **No held-out evaluation ever runs.** The pipeline produces features and stops, so there is no
  number of any kind to discuss.
- **Leakage.** Admixtures sharing a donor or tumor source land in both train and test, so whatever
  number is reported is meaningless. The admixture plan's orthogonality constraint between source
  membership and `f` exists for this reason; the split must honour it.
- **No ichorCNA head-to-head.** The entire premise is *complementary to*. Without the comparison on
  the same samples there is no story regardless of the model's own numbers.
- **Stalling mid-infrastructure.** The honest risk for a project of this shape is dying at "multi-TB
  BWA-MEM on a shared cluster" rather than at "the model didn't work." A Snakemake pipeline with no
  evaluation attached is not a result.

## 4. What this implies for scope

Scope toward **any defensible number, fast** — then improve it. Prefer a deliberately crude
end-to-end result on a small subset (HCC1395 + a handful of healthy samples, coarse bins, no GC
correction) over a full-quality pipeline evaluated at the very end. The crude run buys the thing
that makes the project a result at all: a train/test split with a number on it and an ichorCNA
column next to it. Everything after that is refinement, and refinement narrates cleanly.

Concretely, this ordering principle should be applied when triaging
the backlog: work that shortens the path to a first held-out number outranks work
that improves the quality of a number that does not yet exist.

## 5. How the result gets written up

The writeup must state the **hypothesis, the design, and the decision rule before the results
section** — including the tumor-fraction floor and error metric that would count as success, and the
ichorCNA comparison the result is judged against. Fixing those in advance is what makes a negative
result cost nothing: the reader can see the null would have been reported just as readily as a win.