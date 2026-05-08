# Pipeline fixes — end-to-end validation on 2026-05-02 / 2026-05-03

This document logs every issue found while validating the full TAX + MAG
pipeline (and the `--skip_gtdbtk` Sourmash variant) on a fresh checkout
of the `dev` branch, plus the fix applied to each one.

The intent is to give a self-contained record that can be folded into
the eventual commit message or kept as an internal post-mortem.

## Final status

Both variants now run end-to-end with `exit status 0` and produce all
expected outputs:

```bash
# GTDB-Tk variant (full r232 reference, ~200 GB RAM):
nextflow run main.nf \
    --input  /home/asanzc/CENSALUD/DATA/260226_EPIM232 \
    --outdir /home/asanzc/CENSALUD/RESULTS/260226_EPIM232 \
    --db_root /home/asanzc/CENSALUD/DATABASES \
    --run_assembly \
    -profile gru -resume
# → 32_reports/<run>_epitaxmag_report.html (≈ 5 MB)
# → 32_reports/<run>_full_report.xlsx
# → 27_amr_pathogen_integration/<sample>_{amr_pathogen.tsv,summary.txt}

# Sourmash variant (lightweight, <8 GB RAM):
nextflow run main.nf \
    --input  /home/asanzc/CENSALUD/DATA/260226_EPIM232 \
    --outdir /home/asanzc/CENSALUD/RESULTS/260226_EPIM232 \
    --db_root /home/asanzc/CENSALUD/DATABASES \
    --run_assembly --skip_gtdbtk \
    -profile gru -resume
# Same outputs, plus 24_taxonomy_sourmash/<sample>/
```

## Fixes index

| # | Area                       | What broke                                      |
|---|----------------------------|-------------------------------------------------|
| 1 | nextflow.config            | Per-tool `conda` lost when profile redefines `withName:` |
| 2 | workflow/mag/main.nf       | Final reports defined but never invoked          |
| 3 | main.nf                    | `combine(ch_setup_ready)` invalidated `-resume` cache |
| 4 | scripts + envs/nf-html-report.yml | `os.makedirs("")` + missing `openpyxl`     |
| 5 | workflow/mag/main.nf       | CHECKM2 join key was the work-dir hash, not sample id |
| 6 | workflow/mag/main.nf       | MAG_REPORT cached-on-constant trigger never re-ran |
| 7 | nextflow.config + workflow | Skani as the new default MAG taxonomy classifier  |
| 8 | workflow/mag/main.nf       | Sourmash gather-headers concatenation bug         |
| 9 | conf/{low,medium,high}     | Replace gru/pc/censalud profiles by size-tier ones |
|10 | conf/<size> + main.nf      | Tighten per-process CPU/RAM and raise `maxForks`  |
|11 | workflow/mag/main.nf       | AMRFinderPlus DB path was one level too high     |
|12 | nextflow.config            | Bakta DB path missed the inner versioned directory |
|13 | scripts/*                  | Report and integrator hardcoded to GTDB-Tk format |
|14 | workflow/mag/main.nf + scripts | New AMR_ON_PLASMIDS process (independent AMRFinderPlus pass) |



## Fix 1 — `withName:` blocks inside profiles wipe out conda directives

### Symptom

`EPITAXMAG_MAG:METAFLYE` exited with status 127 and stderr:

```
.command.sh: line 2: flye: command not found
```

The conda env activated by `.command.run` was the global default
`nf-metatax.yml`, not `nf-mag-flye.yml`.

### Root cause

`nextflow.config` declared per-tool conda environments at the **global
process scope**:

```groovy
process {
    conda = "${projectDir}/envs/nf-metatax.yml"
    withName: 'METAFLYE' { conda = "${projectDir}/envs/nf-mag-flye.yml" }
    withName: 'MEDAKA'   { conda = "${projectDir}/envs/nf-mag-medaka.yml" }
    // ... etc
}
```

Every hardware profile (`gru`, `pc`, `censalud`) then redeclared the
same selectors to set `cpus`, `memory`, `maxForks`:

```groovy
profiles {
    gru {
        process {
            withName: 'METAFLYE' { cpus = 48; memory = '80 GB'; maxForks = 1 }
        }
    }
}
```

When two `withName:` blocks target the same selector at different
scopes, Nextflow does **not** merge their directives — the more
specific (profile) scope replaces the entire block. `conda` is
therefore lost as soon as a hardware profile is chosen, and the
process falls back to the global default `nf-metatax.yml`.

`nextflow config -profile gru` confirmed this: `withName:METAFLYE`
showed only `cpus`/`memory`/`maxForks`, no `conda`.

### Fix

Move every per-tool `conda = ...` directive **into each profile's
`withName:` block**, alongside `cpus`/`memory`. The global `process`
scope keeps only the default `conda = nf-metatax.yml` (used by TAX
processes that don't override it) plus any selectors that don't have
profile-specific resource overrides (e.g. `FINAL_REPORT`,
`AMR_PATHOGEN_INTEGRATION|MAG_REPORT`).

This duplicates the conda paths across the three hardware profiles,
but it is the smallest change that respects Nextflow's actual scope
resolution rules and keeps the file readable.

After the fix, `nextflow config -profile gru` shows
`conda = ...nf-mag-flye.yml` for `withName:METAFLYE`, and the process
activates the right env at runtime.


## Fix 2 — final reports never produced because workflow MAG never invokes them

### Symptom

The pipeline ran end-to-end with exit status 0, but the output directory
`32_reports/` was never created and the integration TSV
`27_amr_pathogen_integration/` was missing too. The Nextflow log never
showed `Submitted process > EPITAXMAG_MAG:AMR_PATHOGEN_INTEGRATION` or
`MAG_REPORT`.

### Root cause

`workflow/mag/main.nf` defines both `AMR_PATHOGEN_INTEGRATION` and
`MAG_REPORT` as processes but the `workflow MAG { … }` body never calls
them. The DAG terminates after `BAKTA`/`INTEGRONFINDER`, so Nextflow
considers the workflow complete and exits cleanly.

The downstream processes also could not be wired up trivially: the
upstream emitters (`GTDBTK`, `SOURMASH_CLASSIFY`, `AMRFINDERPLUS_MAGS`,
`CHECKM2`, `GENOMAD`) emit bare `path(…)` channels rather than
`tuple(val(sample), path(…))`, so a multi-source `join` to feed
`AMR_PATHOGEN_INTEGRATION` (one task per sample) is awkward.

### Fix

Two changes in `workflow/mag/main.nf`:

1. Make the upstream MAG processes emit `tuple val(sample), path(…)`
   for the artefacts that downstream tasks consume. This is a no-op
   for caching (Nextflow keys on inputs, not output shape) and keeps
   each artefact paired with its sample for `join`.
   Affected processes:
   - `CHECKM2` (`report`)
   - `GTDBTK` (`summary`)
   - `SOURMASH_CLASSIFY` (`summary`)
   - `AMRFINDERPLUS_MAGS` (`report`)
   - `GENOMAD` (`plasmid_summary`)

2. Append two steps to the `workflow MAG { … }` body:
   - `AMR_PATHOGEN_INTEGRATION`: built by `join`-ing the four per-sample
     channels above (taxonomy, AMR, CheckM2, plasmids). Uses
     `SOURMASH_CLASSIFY.out.summary` when `--skip_gtdbtk` is set.
   - `MAG_REPORT`: a single trigger task that fires once everything
     else has emitted, by mixing all per-sample outputs into a single
     `ready` channel.

After the fix, the pipeline produces `27_amr_pathogen_integration/` and
`32_reports/<run>_epitaxmag_report.html` + `<run>_full_report.xlsx`.

## Fix 3 — `combine(ch_setup_ready)` gating breaks `-resume` caches

### Symptom

After successfully completing a TAX+MAG run, re-launching with
`-resume` re-ran every TAX task from scratch (PORECHOP, CHOPPER,
KRAKEN2, KAIJU, METAFLYE, ...). The cache hit ratio collapsed from
~all-cached to "Submitted" everywhere despite the work directories
existing on disk.

### Root cause

`main.nf` previously inserted a synthetic dependency between
`SETUP_DB` and the downstream subworkflows by mixing every setup
output into a single `ready` signal and combining it with the input
channel:

```groovy
ch_tax_input = ch_raw_fastq
    .combine(ch_setup_ready)
    .map { sample, fastq, _ready -> tuple(sample, fastq) }
EPITAXMAG_TAX(ch_tax_input)
```

The intent was sound (don't start TAX before its DBs exist), but the
side-effect is that Nextflow records `ch_setup_ready` as part of the
input fingerprint of every TAX task. Whenever any setup process had a
new output declaration or work-dir hash (which happened naturally as
we declared more files as outputs of `DOWNLOAD_KAIJU` /
`DOWNLOAD_KRAKEN2`), the ready signal's hash changed and **every
single downstream task lost its cache**, including the expensive ones
like METAFLYE, MEDAKA and GTDB-Tk.

### Fix

Drop the explicit `combine(...)` gating from `main.nf`. The setup
processes already use `storeDir`, so they are idempotent — a
re-launch is a near-instant cache-hit when the DBs are already on
disk.

For the genuine first-run race (TAX:FASTQSCREEN can start before
SETUP_FASTQSCREEN finishes downloading the 1.1 GB Zenodo panel) we
rely on Nextflow's `errorStrategy = 'retry'` (default for the
pipeline) plus the documented operator workflow: if the first launch
trips on a missing DB file, simply relaunch with `-resume` once the
storeDir is populated. This is documented in the `main.nf` workflow
comment and in the README quickstart.

The trade-off (longer first-run setup window with a possible retry)
is much smaller than the cost of invalidating the resume cache for
every downstream MAG task on every config change.


## Fix 4 — MAG_REPORT crash: empty dirname + missing openpyxl

### Symptoms

`MAG_REPORT` failed at the very end of the pipeline with:

```
FileNotFoundError: [Errno 2] No such file or directory: ''
  ... os.makedirs(os.path.dirname(args.output), exist_ok=True)
ModuleNotFoundError: No module named 'openpyxl'
```

### Root cause

1. `scripts/generate_mag_report_v2.py` does
   `os.makedirs(os.path.dirname(args.output), exist_ok=True)`. When the
   process runs in a Nextflow work dir, `args.output` is a bare
   filename like `<run>_epitaxmag_report.html` and `os.path.dirname()`
   returns the empty string, which `os.makedirs("")` rejects.
2. `envs/nf-html-report.yml` ships pandas + plotly but not `openpyxl`,
   which `scripts/generate_excel_report.py` imports at module load.

### Fix

1. In `generate_mag_report_v2.py`, replace
   `os.makedirs(os.path.dirname(args.output), exist_ok=True)` with
   `os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)`.
   `generate_excel_report.py` already had this guard.
2. Add `openpyxl>=3.1` to `envs/nf-html-report.yml`. Nextflow detects
   the env-yml change and rebuilds the conda env on next run.


## Fix 5 — AMR_PATHOGEN_INTEGRATION never fires because CHECKM2's join key is wrong

### Symptom

After Fix 2 wired AMR_PATHOGEN_INTEGRATION into the workflow, the
process still never showed up in the Nextflow dispatch log and
`27_amr_pathogen_integration/` was never produced. The join silently
emitted nothing.

### Root cause

The four channels feeding `AMR_PATHOGEN_INTEGRATION` had to be tagged
with the sample id. Three of them (taxonomy, AMR, plasmids) carry the
sample id in the filename and could be derived with a regex on
`f.baseName`. CHECKM2 is the exception — it emits a generic
`quality_report.tsv` (the sample is only encoded in the `publishDir`
path, not in the filename or the channel item).

I had reached for `.map { f -> tuple(f.parent.name, f) }`, expecting
`f.parent.name` to be the publish directory's sample folder. It is
not: `f.parent` resolves to the Nextflow work directory (e.g.
`work/24/0cd729…`) and `.name` is the random hash. The `join` then
matched a hash key against the real sample id and silently produced
zero items.

### Fix

Make `CHECKM2` emit `tuple val(sample), path("quality_report.tsv")`
instead of bare `path(...)`. This is the same kind of tuple shape we
deliberately avoided for the bigger upstream emitters (to keep their
caches valid) — but for CHECKM2 the cost is a single ~5-minute
re-run, and downstream `AMR_PATHOGEN_INTEGRATION` then joins cleanly
because every channel uses the same key (the sample id) without any
filename parsing.


## Fix 6 — MAG_REPORT served stale HTML when AMR_PATHOGEN_INTEGRATION ran later

### Symptom

After Fix 5 made `AMR_PATHOGEN_INTEGRATION` finally fire on a `-resume`
run, the integration TSV was produced under
`27_amr_pathogen_integration/` but the final
`32_reports/<run>_epitaxmag_report.html` was not regenerated; it was
served from the previous run's cache and therefore did not include
the AMR-pathogen sections.

### Root cause

`MAG_REPORT` was wired with `input: val(ready)` and triggered with a
synthetic constant string `'ready'`. Nextflow hashes inputs to detect
work that needs redoing — and the constant `'ready'` always hashes to
the same value. So as soon as `MAG_REPORT` had been executed once it
was a permanent cache hit, regardless of whether any of the artefacts
the report reads from `${outdir}` had since changed.

### Fix

Change `MAG_REPORT`'s input to `path(triggers)` and pass the actual
`AMR_PATHOGEN_INTEGRATION.out.report`, `BAKTA.out.all_files`,
`INTEGRONFINDER.out.summary`, `ABRICATE.out.vfdb` and
`MOBSUITE.out.results` files (mixed and collected) as the upstream
channel. Nextflow now hashes those staged files, so any change to
them invalidates the cache and the report is rebuilt. The script
itself still reads from `${outdir}` — the staged files act purely as
a synchronization barrier and a hash signal.


## Fix 7 — Skani is now the default MAG taxonomy classifier

### Motivation

GTDB-Tk r232 is the gold-standard but needs ≥120 GB RAM (pplacer is
the bottleneck) and a ~60 GB reference package. That made the default
pipeline unusable on the deployment target (≤64 GB RAM, El Salvador
collaborators) and turned every `-resume` cycle that touched
`workflow/mag/main.nf` into a multi-hour GTDB-Tk re-run.

Skani v0.3 produces real ANI vs the GTDB r226 sketch DB, fits in
≤24 GB RAM, runs in seconds per MAG, and the bluenote-1577 lab ships
a pre-built single-file sketch (`skani_gtdb_r226-v0.3.tar.gz`,
~38 GB compressed → ~36 GB extracted) ready to plug in. For
surveillance use-cases (water, AMR genus/species ID) Skani's accuracy
is on par with GTDB-Tk.

### Changes

- New `params.taxonomy_tool` selector with three values: `skani`
  (default), `sourmash`, `gtdbtk`. Replaces the boolean
  `params.skip_gtdbtk`, which is kept as a deprecated alias that
  rewrites the selector to `sourmash` (with a warning) so existing
  scripts keep working.
- New `envs/nf-mag-skani.yml` (skani ≥0.3 + python).
- New `DOWNLOAD_SKANI` process in `workflow/setup_db.nf` — fetches
  the pre-built sketch into `${db_root}/skani/skani_gtdb_r226-v0.3/`.
- New `SKANI_CLASSIFY` process in `workflow/mag/main.nf` — runs
  `skani search -d <db> -q <bins>` once per sample, then converts
  the result into a GTDB-Tk-shaped TSV (one best-ANI row per bin)
  so that `AMR_PATHOGEN_INTEGRATION` does not need to know which
  classifier produced it.
- `workflow MAG` and `workflow SETUP_DB` switch on
  `params.taxonomy_tool` to wire Skani / Sourmash / GTDB-Tk.
- README and CLAUDE.md updated: the RAM table now shows Skani as
  the default, GTDB-Tk demoted to opt-in.

### Trade-offs

- Pipeline default RAM bar drops from ~200 GB to ≤24 GB.
- DB footprint drops from ~60 GB (GTDB-Tk r232) to ~36 GB (Skani r226).
- Resolution stops one GTDB release short (r226 vs r232) — irrelevant
  for surveillance of canonical waterborne pathogens, but worth noting
  for publications that demand the latest reference.


## Fix 8 — Sourmash classification silently produced an empty taxonomy TSV

### Symptom

When the pipeline ran with the Sourmash branch, `sourmash gather`
emitted hits for every MAG (good), but the per-sample
`<sample>_taxonomy.tsv` came back empty (header only). MAG_REPORT then
showed every SemiBin "without taxonomy".

### Root cause

The post-processing step concatenated the per-bin gather CSVs with:

```bash
head -1 gather_results/*.gather.csv | head -1 > gather_all.csv
```

When `head -1` is given multiple files it prefixes each with a
`==> file <==` banner, so the first line piped into the inner
`head -1` is a banner, not the real CSV header. `sourmash tax genome`
then refused the malformed CSV and the script silently fell through
(the `2>/dev/null || true` swallowed the error), leaving the output
TSV empty.

### Fix

Pick a single CSV explicitly to extract the header instead of relying
on `head -1` over a glob:

```bash
first_gather=$(ls gather_results/*.gather.csv 2>/dev/null | head -1)
if [ -n "$first_gather" ]; then
    head -1 "$first_gather" > gather_all.csv
    for f in gather_results/*.gather.csv; do
        tail -n +2 "$f" >> gather_all.csv
    done
fi
```


## Fix 9 — Replace machine-specific profiles with size-tier profiles

### Motivation

The original `gru` / `pc` / `censalud` profiles named the deployment
hosts the developer happened to use, which (a) leaks internal hostnames
into the public README and (b) tied future deployments to whatever
hardware those names implied. New collaborators have to read a wiki
to know whether `pc` or `gru` matches their workstation.

### Changes

- New `conf/low.config`, `conf/medium.config`, `conf/high.config`
  with explicit hardware targets:
  * `low`     ≥16 CPU / ≥32 GB RAM
  * `medium`  ≥24 CPU / ≥64 GB RAM (recommended default)
  * `high`    ≥64 CPU / ≥128 GB RAM (only profile that allows
    `--taxonomy_tool gtdbtk`).
- The `nextflow.config` `profiles { … }` block shrinks to a thin
  registry that `includeConfig`s the right size file. `pc`,
  `censalud` and `gru` are kept as aliases that include the matching
  size config; `main.nf` logs a deprecation warning when one of the
  legacy names is selected.
- `main.nf` rejects `--taxonomy_tool gtdbtk` on any profile other
  than `high` / `gru`, so users cannot accidentally OOM their host
  by asking pplacer to fit in 64 GB.


## Fix 10 — Tighten per-process CPU / RAM and raise `maxForks`

### Motivation

The previous resource blocks were generous to the point of waste —
KRAKEN2 booked 48 CPUs and 120 GB on every host, BAKTA asked for
16 CPUs even when annotating a 4 MB MAG, and `maxForks = 1` was the
default for every binner. With a single sample that is harmless;
with a typical 7-sample Nanopore run, the binning stage was running
each binner sequentially (≈ 14 h) when it could have been running
two or three at once (≈ 5 h).

### Changes per category (numbers given for `medium`)

- **Light Python / Java** (NANOPLOT, FASTQC, MULTIQC, FINAL_REPORT,
  MAG_REPORT, AMR_PATHOGEN_INTEGRATION, BRACKEN, PORECHOP, CHOPPER,
  CHECK_SQCORE_*): 2–4 CPUs, 4–8 GB; `maxForks` left unbounded.
- **Mapping family** (FASTQSCREEN, HOST_REMOVAL, MAP_READS,
  KMA_READS, MINIMAP2_PHENOTYPIC): 12 CPUs, 16 GB; `maxForks = 3`.
- **Heavy k-mer classifiers** (KRAKEN2, KAIJU): 20 CPUs, 20 GB,
  always memory-mapped on `low`/`medium`; `maxForks = 1`.
- **Assembly** METAFLYE: 20 CPUs, 48 GB; `maxForks = 1`.
- **Polishing** MEDAKA: 2 CPUs, 6 GB; `maxForks = 8` (single-thread
  NN, parallelizing over samples is the only way to fill the host).
- **Binners** METABAT2 / MAXBIN2 / SEMIBIN2: 8–12 CPUs, 16–20 GB;
  `maxForks = 2`.
- **MAG QC + taxonomy** CHECKM2 12 CPUs / 20 GB / `maxForks = 2`,
  SKANI_CLASSIFY 12 CPUs / 24 GB / `maxForks = 2`,
  SOURMASH_CLASSIFY 6 CPUs / 8 GB / `maxForks = 6`.
- **Per-bin annotation** BAKTA, INTEGRONFINDER, AMRFINDERPLUS_MAGS,
  ABRICATE: 6 CPUs, 6–8 GB; `maxForks = 4–6` so multiple bins of the
  same sample run in parallel inside a single sample.

### Safety net

The per-process `maxForks` is always paired with the executor caps
(`executor.cpus` / `executor.memory` per profile). Even if every
binner could in theory launch its `maxForks` quota, Nextflow refuses
to schedule a task once the running tasks would exceed the executor's
total CPU or memory budget. This is what prevents two MetaFlye + a
SemiBin2 pile-up on the same host.


## Fix 11 — AMRFinderPlus per-MAG silently produced empty output

### Symptom

After several end-to-end runs, the AMR sections of the report (MAG
table, AMR-by-class heatmap, MAG ↔ AMR cross-reference) were always
empty. KMA-on-reads found 19 hits and geNomad reported plasmids with
ARGs, so the data clearly existed somewhere — yet the per-MAG TSV
`25_amr_mags/<sample>/<sample>_amr_combined.tsv` and the directory
`per_bin/` came back empty.

### Root cause

Two compounding bugs:

1. The DB path passed to `amrfinder --database` pointed at
   `${db_root}/amrfinderplus/latest`. `amrfinder_update` lays out the
   data as `<db_root>/amrfinderplus/latest/<version>/<files>` plus a
   convenience symlink `<db_root>/amrfinderplus/latest/latest →
   <version>/`. The pipeline was passing the parent of that symlink,
   so `amrfinder` looked for `AMRProt.fa.phr` directly inside `latest/`
   and bailed with `*** ERROR *** The BLAST database for AMRProt.fa
   was not found.`
2. Stderr was redirected to `/dev/null` and the call ended in
   `|| true`, so the BLAST-DB-not-found message never reached the
   Nextflow log. Every per-MAG invocation failed silently and produced
   no `.tsv`, leaving an empty directory.

### Fix

- Add `params.amrfinder_db = "${db_root}/amrfinderplus/latest/latest"`
  in `nextflow.config` (with a comment about the symlink layout) and
  switch the `AMRFINDERPLUS_MAGS` script to use it.
- Redirect stderr to a per-bin log file (`per_bin/<bin>.log`) instead
  of `/dev/null` so future DB / index issues surface immediately.
- Replace the `head -1 *.tsv` concatenation with the explicit
  single-file pattern already used elsewhere (avoids the same
  `==> file <==` banner pitfall as Sourmash had).

After the fix, AMRFinderPlus locates the formatted BLAST indices and
finds, for example, a `blaOXA-211` carbapenem-hydrolyzing class D
β-lactamase in the *Acinetobacter johnsonii* MAG — the kind of hit
the report was meant to surface.


## Fix 12 — Bakta DB path missed the inner versioned directory

### Symptom

`Bakta summaries: 0` in MAG_REPORT regardless of run. Per-MAG `.bakta.log`
contained:

```
ERROR: database version file (version.json) not readable!
Please check if /home/asanzc/CENSALUD/DATABASES/bakta/db is the correct
path. Maybe there is another 'db' subdirectory?
```

### Root cause

Bakta's installer extracts the database into `<root>/bakta/db/db/`
(`version.json`, `bakta.db`, `amrfinderplus-db/`, etc. all live in the
inner `db/`). `params.bakta_db` was pointing one level too high.

### Fix

Set `bakta_db = "${params.db_root}/bakta/db/db"` in `nextflow.config`
and document the layout inline.


## Fix 13 — Report assumed GTDB-Tk taxonomy and integrator could not parse Skani

### Symptom

After switching the default to Skani (Fix 7), the HTML report still
said “GTDB-Tk r232” everywhere, and `AMR_PATHOGEN_INTEGRATION` filled
the `Organism` column with `Unknown` because `parse_taxonomy_levels`
only understood the `g__/s__` GTDB notation. The MAG catalogue also
showed no Bakta-derived stats (no CDS / tRNA / rRNA / CRISPR columns)
and the geNomad `amr_genes` column was hidden behind a per-sample
summary that collapsed it.

### Fix

`scripts/generate_mag_report_v2.py`:

- New `detect_taxonomy_tool(rd)` that resolves the `24_taxonomy_*`
  directory actually populated and returns label / version / db tag.
  All hardcoded "GTDB-Tk r232" strings (DAG box, software list, MAG
  catalogue caption, taxonomy-section caption, log line) now read
  from this struct, so swapping `--taxonomy_tool` reflects in the
  report without further edits.
- `load_taxonomy()` parses both formats: GTDB rank prefixes when the
  directory is `24_taxonomy_gtdbtk`, and a free-text best-hit string
  (`<accession> <Genus species ...>`) for Skani / Sourmash. The
  "Class" column used by the per-sample stacked bar falls back to
  genus when no GTDB class rank is available.
- MAG catalogue table now exposes Bakta CDS, tRNA, rRNA, ncRNA,
  CRISPR counts and a per-MAG AMR-gene count.
- Plasmid section reorganised into a per-sample summary plus a
  per-contig detail of confirmed plasmids carrying AMR (taken from
  geNomad's built-in `amr_genes` column).
- New “All KMA/ResFinder hits in reads” table inside the concordance
  section, plus a placeholder for the AMR-on-plasmids panel filled
  in by Fix 14.

`scripts/integrate_amr_pathogen.py`:

- `parse_taxonomy_levels()` accepts both formats: detects GTDB by the
  presence of `d__/p__/g__/s__` tokens, otherwise treats the string as
  a Skani / Sourmash best hit and extracts genus / species from the
  binomial in the title up to common stop-tokens (`strain`, `isolate`,
  `contig`, `,` …).


## Fix 14 — AMRFinderPlus on geNomad plasmid contigs (independent pass)

### Symptom

The plasmid section showed only the AMR genes that geNomad annotated
itself (using its bundled MMseqs2 profile DB). That captures core ARGs
but misses anything outside the geNomad model. Worse, plasmid contigs
that never made it into a high-quality MAG were invisible to the
AMR-in-MAGs pipeline.

### Fix

New `AMR_ON_PLASMIDS` process in `workflow/mag/main.nf`:

- `GENOMAD` now also emits `${sample}_plasmid.fna` (was internal-only).
- `AMR_ON_PLASMIDS(GENOMAD.out.plasmid_fna)` runs `amrfinder
  --nucleotide` on the plasmid FASTA, reusing
  `params.amrfinder_db` from Fix 11.
- Output is reshaped on the fly into the columns the report expects
  (`Sample / Contig / Gene / Class / Subclass / Identity / Coverage /
  Accession / Closest_ref`) and published to
  `26_genomad/<sample>/<sample>_plasmid_amr.tsv`.
- An empty plasmid FASTA produces a header-only TSV so downstream
  steps stay happy.

`scripts/generate_mag_report_v2.py` adds `load_amr_plasmids()` and a
new sub-section "AMRFinderPlus on plasmid sequences" inside section 13
that renders the hits with their identity / coverage / class /
subclass / NCBI-protein link.
