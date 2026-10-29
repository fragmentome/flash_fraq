# Per-workflow Snakemake working directories

Snakemake keeps its state -- the directory lock, job metadata, incomplete
markers, and the SLURM executor's per-job logs -- in a `.snakemake`
directory inside whatever working directory it was invoked from. With every
driver invoking it from the repo root, all of MAT-19/20/23 shared one
`.snakemake`, so a second workflow could not start while a first was
running: it died at DAG build on the directory lock. MAT-20 was blocked this
way behind the 21-day MAT-23 cohort run (job 2749087).

Each cohort/smoke driver therefore runs Snakemake from its own subdirectory
here, giving it its own `.snakemake` while still using the single repo
checkout for code.

    snakemake_workdirs/
      mat19/          .snakemake for scripts/drivers/mat19_driver.sh
      mat20/          ... mat20_driver.sh
      mat23/          ... mat23_driver.sh
      hcc1395_ena_local/
                      ... hcc1395_ena_local_driver.sh (workstation, one run
                          from manifests/hcc1395/ straight off ENA)
      ...             one per driver, named after its RUN_NAME

The subdirectories are created by the drivers on demand and their contents
are gitignored -- only this README is tracked.

Two things this requires, both already in place:

* **Nothing in the workflow may use a path relative to the invocation
  directory.** Snakemake resolves a relative `configfile:` or `shell:` path
  against the working directory, so `python scripts/foo.py` would look for
  the script inside `snakemake_workdirs/<name>/`. `workflow/Snakefile`
  derives `REPO_ROOT` from `workflow.basedir` and every rule builds its
  script path from `SCRIPTS_DIR`.
* **`-s` and `--workflow-profile` must be absolute** on the command line,
  since they are resolved before Snakemake changes directory.

A consequence worth knowing: state is per directory, so `--unlock`,
`--rerun-incomplete` and `--cleanup-metadata` only ever affect the workflow
whose directory you run them in. `scripts/drivers/_preflight.sh` relies on
this -- its gate now only blocks on jobs belonging to the same workflow,
because the drivers submit with `sbatch --chdir=<this workflow's dir>`, so
both the orchestrator and every per-rule job the executor plugin submits
record that directory as their SLURM WorkDir.
