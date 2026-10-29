# flash fraq documentation index

This index covers everything under `docs/`.

## Project design and plan

| Doc | What's in it |
| --- | --- |
| [motivation-and-goal.md](motivation-and-goal.md) | Why the project exists and what counts as finishing it — the falsifiable question, why a negative result is a complete result, the four failure modes that would actually sink it, and the scope/writeup rules that follow |
| [phase0-spec.md](phase0-spec.md) | The technical spec: ichorCNA vs. DELFI grounding, the structural gap this project targets, feature set, in-silico-admixture labeling strategy |
| [phase1-sourcing.md](phase1-sourcing.md) | Data sourcing verification — which cohorts exist, what's actually reachable, how to get SEQC2 (HCC1395/HCC1395BL). Some "already downloaded" claims here proved stale; see the TODO files |
| [phase2-preprocessing-plan.md](phase2-preprocessing-plan.md) | End-to-end preprocessing plan: reference build, per-source alignment paths, dedup strategy, fragment-size extraction, GC-bias correction and binning |

## Implementation reference

| Doc | What's in it |
| --- | --- |
| [pipeline-modules.md](pipeline-modules.md) | `src/flash_fraq` module reference and callstack, in data-dependency order — the authority on pipeline order and the "dedup once, post-merge" invariant |
| [local-features-pipeline.md](local-features-pipeline.md) | The workstation feature-extraction DAG: stage-in/analyse/purge, `temp()` as the purge mechanism, the `ssd_mb` concurrency budget |
| [admixture-design.md](admixture-design.md) | The mixture plan: harvested native coverage for all 547 samples, why `C_target` is 0.9x, the fraction ladder, and how the plan keeps source membership orthogonal to `f`. **⚠️ Background pool superseded 2026-08-31** — background is now HCC1395BL, the matched normal; plan must be regenerated |
| [cohort-qc-baseline.md](cohort-qc-baseline.md) | First run of the pre-admixture cohort matrix + phase2 §6 decoy-classifier gate over the real Snyder + Cristiano samples, and why its `passed: true` is necessary but not sufficient |

## Infrastructure

| Doc | What's in it |
| --- | --- |
| [workstation-ssd-scratch.md](workstation-ssd-scratch.md) | Why CRAMs are staged on SSD before analysis and how the size constraint is enforced |

