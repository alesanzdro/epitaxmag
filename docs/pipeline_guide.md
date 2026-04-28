# EpiTaxMAG — Pipeline Documentation

**Version:** 1.0  
**Date:** April 2026  
**Authors:** Alejandro Sanz Carbonell — FISABIO EPIMOL (Molecular Epidemiology), Valencia, Spain  
**Application:** Metagenomic surveillance of water samples using Oxford Nanopore sequencing  
**License:** GPL-3.0

---

## Table of Contents

1. [Architecture](#1-architecture)
2. [Phase TAX (EpiTax) — Steps 01–14](#2-phase-tax-epitax--steps-0114)
3. [Phase MAG (EpiTaxMAG) — Steps 15–32](#3-phase-mag-epitaxmag--steps-1532)
4. [Phenotypic Verification Module](#4-phenotypic-verification-module)
5. [Reports](#5-reports)
6. [Databases](#6-databases)
7. [Hardware Profiles](#7-hardware-profiles)
8. [Singularity & Conda](#8-singularity--conda)
9. [Key Design Decisions](#9-key-design-decisions)
10. [Process Reference Table](#10-process-reference-table)

---

## 1. Architecture

EpiTaxMAG is a modular Nextflow DSL2 pipeline with two independent but complementary phases:

- **EpiTax** (TAX phase, steps 01–14): Quality control, adapter trimming, quality filtering, optional host depletion, multi-tool taxonomic profiling, read-level AMR screening, phenotypic verification, and interactive HTML report.

- **EpiTaxMAG** (MAG phase, steps 15–32): Metagenomic assembly, multi-algorithm binning, MAG quality assessment, taxonomic classification, functional annotation, viral/plasmid detection, integron analysis, and comprehensive reports with AMR-pathogen-plasmid risk assessment.

### Execution modes

```bash
# TAX only (default) — ~4 hours for 7 samples
nextflow run main.nf --input data/EPIM232 -profile gru

# TAX + MAG — ~12-24 hours
nextflow run main.nf --input data/EPIM232 -profile gru --run_assembly

# MAG only (reuses filtered reads from a prior TAX run)
nextflow run main.nf --input data/EPIM232 -profile gru --run_tax false --run_assembly
```

### Nextflow process naming

All processes appear with the `EPITAXMAG_` prefix in the Nextflow log, following nf-core convention:

```
EPITAXMAG_TAX:KRAKEN2        (EPI001155_PA004)  [100%] 7 of 7
EPITAXMAG_TAX:BRACKEN        (EPI005151_100_2e) [100%] 7 of 7
EPITAXMAG_MAG:METAFLYE       (EPI001155_PA004)  [100%] 7 of 7
```

This is achieved via aliased workflow imports in `main.nf`:

```groovy
include { TAX as EPITAXMAG_TAX } from './workflow/tax/main'
include { MAG as EPITAXMAG_MAG } from './workflow/mag/main'
```

---

## 2. Phase TAX (EpiTax) — Steps 01–14

All processes defined in `workflow/tax/main.nf`.

### Step 01 — QC Raw Reads
**Directory:** `01_qc_raw/`  
**Processes:** `NANOPLOT_RAW`, `FASTQC_RAW`, `CHECK_SQCORE_RAW`  
**Tools:** NanoPlot 1.46.2, FastQC 0.12.1, custom Python (Dorado metrics)  
**Description:** Quality metrics on raw reads: read length distribution, N50, GC content, quality scores. CHECK_SQCORE extracts Dorado-specific metrics (Q-score, basecalling model) into a CSV with European decimal format (`;` separator, `,` decimal).

### Step 02 — FastQ Screen
**Directory:** `02_fastqscreen_raw/`  
**Process:** `FASTQSCREEN`  
**Tools:** FastQ Screen 0.16.0, minimap2 2.30  
**Description:** Contamination screening against a custom 30-genome panel (host, controls, fecal indicators, pathogens, viruses, parasites). Runs on **raw** reads (before trimming) to detect vectors and adapters. Uses `--aligner minimap2` with `--minimap2 '-x map-ont'` for Nanopore compatibility. Panel hosted on Zenodo for automatic download.

### Step 03 — Adapter Trimming
**Directory:** `03_trim_porechop/`  
**Process:** `PORECHOP`  
**Tools:** Porechop ABI 0.5.1  
**Description:** Adapter and chimera removal from Nanopore reads. Porechop ABI is the actively maintained fork.

### Step 04 — Quality & Length Filtering
**Directory:** `04_filter_chopper/`  
**Process:** `CHOPPER`  
**Tools:** Chopper 0.12.0  
**Description:** Filters reads by quality (Q >= 10) and length (>= 1000 bp). These thresholds are configurable via `--min_quality` and `--min_length`.

### Step 05 — Host Removal (optional)
**Directory:** `05_host_removal/`  
**Process:** `HOST_REMOVAL`  
**Tools:** minimap2 2.30, samtools 1.21  
**Description:** Depletes host reads by mapping against a reference genome and extracting unmapped reads (`samtools view -f 4`). Only runs if `--host_genome /path/to/genome.fasta` is provided. For human samples: GRCh38.

### Step 06 — QC Filtered Reads
**Directory:** `06_qc_filtered/`  
**Processes:** `NANOPLOT_FILTERED`, `FASTQC_FILTERED`, `CHECK_SQCORE_FILTERED`  
**Description:** Same QC as step 01 but on filtered reads, enabling before/after comparison.

### Step 07 — Kraken2 Classification
**Directory:** `07_tax_kraken2/`  
**Process:** `KRAKEN2`  
**Tools:** Kraken2 2.17.1  
**Database:** PlusPF (2025-10-15), ~70 GB  
**Description:** k-mer based taxonomic classification. Uses `--memory-mapping` automatically when RAM < 100 GB. Runs with `maxForks=1` due to high RAM usage.

### Step 08 — Bracken Abundance Estimation
**Directory:** `08_tax_bracken/`  
**Process:** `BRACKEN`  
**Tools:** Bracken 3.1  
**Description:** Re-estimates species-level (S) and genus-level (G) abundances from Kraken2 reports. Corrects for sequence length bias in k-mer classification.

### Step 09 — Kaiju Protein Classification
**Directory:** `09_tax_kaiju/`  
**Process:** `KAIJU`  
**Tools:** Kaiju 1.10.1  
**Database:** RefSeq reference proteomes  
**Description:** Protein-level classification. Complements k-mer approaches by using translated protein sequences. Particularly useful for divergent organisms.

### Step 10 — Sylph ANI Profiling
**Directory:** `10_tax_sylph/`  
**Process:** `SYLPH`  
**Tools:** Sylph 0.9.0  
**Databases:** GTDB r226 + RefSeq Fungi + IMG/VR (viral)  
**Description:** Sketch-based ANI profiling. Provides species-level abundance estimates, effective coverage, and adjusted ANI for each detected organism. Key data source for MAG viability estimation and Whittaker diversity curves.

### Step 11 — AMR Screening (KMA/ResFinder)
**Directory:** `11_amr_kma/`  
**Process:** `KMA_READS`  
**Tools:** KMA 1.6.8  
**Database:** ResFinder (DTU)  
**Description:** Direct read-level detection of antibiotic resistance genes using KMA in `-ont` mode (optimized for Nanopore). No assembly required. Reports gene identity, coverage, and depth. Results classified by antibiotic class (beta-lactam, aminoglycoside, macrolide, tetracycline, quinolone, etc.) in the HTML report.

### Step 12 — Phenotypic Verification (minimap2)
**Directory:** `12_phenotypic_minimap2/`  
**Processes:** `BUILD_PHENOTYPIC_DB`, `MINIMAP2_PHENOTYPIC`  
**Tools:** minimap2 2.30, samtools 1.21, Python  
**Description:** Maps filtered reads against NCBI type-strain reference genomes to verify presence/absence of target organisms. Reports **breadth of coverage** at three depth thresholds: ≥1x (detection), ≥10x (confidence), ≥30x (assembly-ready). Target organisms are defined in `assets/phenotypic_targets/target_organisms.tsv` (editable TSV with organism name + NCBI accession). Genomes are downloaded automatically from NCBI on first run. Enabled by default; disable with `--phenotypic_targets false`.

### Step 13 — MultiQC
**Directory:** `13_multiqc/`  
**Process:** `MULTIQC`  
**Tools:** MultiQC 1.33  
**Description:** Aggregated QC report from FastQC, NanoPlot, Kraken2, Bracken, Kaiju, and KMA outputs.

### Step 14 — EpiTax HTML Report
**Directory:** `14_epitax_report/`  
**Process:** `FINAL_REPORT`  
**Tools:** Python (Plotly, pandas, numpy)  
**Description:** Self-contained interactive HTML report (~5 MB) with embedded Plotly.js. Nine sections:

| Section | Content |
|---------|---------|
| 1. Survival Metrics | Read/base retention before/after filtering |
| 2. FastQ Screen | Contamination screening heatmap |
| 3. Integrated Taxonomy | Per-sample Bracken + Kraken2 + Kaiju charts |
| 4. Diversity Profile | Whittaker rank-abundance curves (Sylph GTDB) |
| 5. MAG Viability | Coverage estimation + depth candidates with 30x threshold |
| 6. AMR Screening | Stacked bars by antibiotic class + gene table + depth |
| 7. Phenotypic Verification | Heatmap + table with breadth @1x/10x/30x |
| 8. Rarefaction Curves | Diversity saturation estimation |
| 9. Conclusion | Automated recommendations |

**Features:** Sortable tables (click header), sample filter checkboxes (hide/show samples across all charts and tables), works offline.

---

## 3. Phase MAG (EpiTaxMAG) — Steps 15–32

All processes defined in `workflow/mag/main.nf`. Activated with `--run_assembly`.

### Step 15 — Metagenomic Assembly
**Directory:** `15_assembly_flye/`  
**Process:** `METAFLYE`  
**Tools:** Flye 2.9.6  
**Description:** De novo metagenomic assembly with `--meta --nano-hq`. The `--nano-hq` preset is correct for Dorado SUP/duplex basecalling. Runs with `maxForks=1` (high memory usage).

### Step 16 — Assembly Polishing
**Directory:** `16_polish_medaka/`  
**Process:** `MEDAKA`  
**Tools:** Medaka 2.2.1  
**Model:** `r1041_e82_400bps_sup_v5.2.0`  
**Description:** Neural network polishing of the Flye assembly. Uses `medaka inference` + `medaka sequence` (Medaka 2.x API). Racon is NOT used — contemporary Medaka models are trained on Flye output directly and Racon can degrade quality with R10.4.1 + SUP v5. Single-threaded inference → 2 CPUs per task for maximum parallelism.

### Step 17 — Assembly QC
**Directory:** `17_qc_assembly/`  
**Process:** `QUAST`  
**Tools:** QUAST 5.3.0  
**Description:** Assembly statistics: N50, total length, number of contigs, GC content.

### Step 18 — Coverage Mapping
**Directory:** `18_coverage/`  
**Process:** `MAP_READS`  
**Tools:** minimap2 2.30, samtools 1.21  
**Description:** Maps filtered reads back to the polished assembly for per-contig depth profiles. Required for MetaBAT2 and SemiBin2 binning.

### Step 19 — MetaBAT2 Binning
**Directory:** `19_binning_metabat2/`  
**Process:** `METABAT2`  
**Tools:** MetaBAT2 2.18  
**Description:** Composition + coverage binning. Requires `jgi_summarize_bam_contig_depths` for depth profiles.

### Step 20 — MaxBin2 Binning
**Directory:** `20_binning_maxbin2/`  
**Process:** `MAXBIN2`  
**Tools:** MaxBin2 2.2.7  
**Description:** Expectation-maximization binning. Designed for short reads but still contributes to consensus.

### Step 21 — SemiBin2 Binning
**Directory:** `21_binning_semibin2/`  
**Process:** `SEMIBIN2`  
**Tools:** SemiBin2 2.2.1  
**Description:** Deep learning binning with `--sequencing-type long_read`. Outputs `.fa.gz` files that are decompressed automatically.

### Step 22 — DAS Tool Consensus
**Directory:** `22_binning_dastool/`  
**Process:** `DAS_TOOL`  
**Tools:** DAS Tool 1.1.7  
**Description:** Bin refinement from MetaBAT2 + MaxBin2 + SemiBin2 consensus. Selects the best non-redundant set of bins. Falls back to single-binner copy if only one produces output. MetaBAT2 scaffolds2bin uses `cut -f1,4` to handle its 4-column format.

### Step 23 — MAG Quality (CheckM2)
**Directory:** `23_binqc_checkm2/`  
**Process:** `CHECKM2`  
**Tools:** CheckM2 1.1.0  
**Description:** Completeness and contamination assessment. Quality tiers: HQ (≥90% completeness, ≤5% contamination), MQ (≥50%, ≤10%), LQ (below MQ).

### Step 24 — MAG Taxonomy
**Directory:** `24_taxonomy_gtdbtk/` or `24_taxonomy_sourmash/`  
**Processes:** `GTDBTK` or `SOURMASH_CLASSIFY`  
**Tools:** GTDB-Tk 2.7.0 (GTDB r232) or Sourmash 4.9.4  
**Description:** Taxonomic classification of MAGs. GTDB-Tk requires >50 GB RAM (pplacer). Sourmash is a lightweight alternative (~5 GB RAM) enabled with `--skip_gtdbtk`. Both use GTDB taxonomy.

### Step 25 — AMR Detection (AMRFinderPlus)
**Directory:** `25_amr_mags/`  
**Process:** `AMRFINDERPLUS_MAGS`  
**Tools:** AMRFinderPlus 4.2.7  
**Description:** Per-MAG detection of AMR genes, virulence factors, and stress response genes. Uses `--plus` for extended detection. Organism-specific analysis when taxonomy is available.

### Step 26 — Viral & Plasmid Detection
**Directory:** `26_genomad/`  
**Process:** `GENOMAD`  
**Tools:** geNomad 1.12.0  
**Description:** Identifies plasmids and viral sequences in the assembly. Critical for determining whether AMR genes are chromosomal or mobile.

### Step 27 — AMR-Pathogen-Plasmid Integration
**Directory:** `27_amr_pathogen_integration/`  
**Process:** `AMR_PATHOGEN_INTEGRATION`  
**Tools:** Python (custom)  
**Description:** Crosses GTDB-Tk taxonomy + AMRFinderPlus + geNomad + CheckM2 to produce AMR-pathogen-plasmid linkage with risk levels: CRITICAL (plasmid-borne AMR in pathogen), HIGH (chromosomal AMR in pathogen), MEDIUM (AMR in unknown host), LOW (AMR in non-pathogen).

### Step 28 — ABRicate
**Directory:** `28_abricate/`  
**Process:** `ABRICATE`  
**Tools:** ABRicate 1.4.0  
**Databases:** VFDB, CARD, PlasmidFinder  
**Description:** Screening against multiple curated databases for virulence factors, resistance genes, and plasmid replicons.

### Step 29 — MOB-suite
**Directory:** `29_mobsuite/`  
**Process:** `MOBSUITE`  
**Tools:** MOB-suite 3.1.9  
**Description:** Plasmid typing, mobility prediction, and replicon classification.

### Step 30 — IntegronFinder
**Directory:** `30_integronfinder/`  
**Process:** `INTEGRONFINDER`  
**Tools:** IntegronFinder 2.0  
**Description:** Detection of class 1/2/3 integrons. Runs on MAGs only (not full assemblies — too slow on 12K+ contigs).

### Step 31 — Bakta Annotation
**Directory:** `31_bakta/`  
**Process:** `BAKTA`  
**Tools:** Bakta 1.12.0  
**Database:** Full (84 GB) or Light (4 GB), configurable via `--bakta_db_type`  
**Description:** Structural and functional annotation: CDS, tRNA, rRNA, CRISPR, signal peptides. Uses `--skip-plot --force` with custom `$TMPDIR` and `--home` for Singularity compatibility.

### Step 32 — MAG Reports
**Directory:** `32_reports/`  
**Process:** `MAG_REPORT`  
**Tools:** Python (Plotly, openpyxl, ReportLab)  
**Description:** Comprehensive reports: interactive HTML (16 sections, sortable/filterable tables, sample filter checkboxes), Excel (16 sheets with all raw data), PDF (formal report with branding).

---

## 4. Phenotypic Verification Module

### Purpose

Verifies the presence of target organisms by mapping filtered reads against NCBI type-strain reference genomes (complete chromosomes + plasmids) using minimap2. Provides quantitative coverage metrics at multiple depth thresholds, enabling confident calls on organism detection.

### Configuration

Target organisms are defined in `assets/phenotypic_targets/target_organisms.tsv`:

```tsv
organism	accession
Escherichia coli	GCF_000005845.2
Klebsiella pneumoniae	GCF_000240185.1
Pseudomonas aeruginosa	GCF_000006765.1
...
```

To customize: edit the TSV — add/remove rows with organism name + NCBI RefSeq accession. Genomes are downloaded automatically from the NCBI Datasets API on first run and cached by Nextflow (`-resume`).

### Metrics reported

| Metric | Meaning |
|--------|---------|
| **Breadth ≥1x** | % of reference genome covered by at least 1 read (detection) |
| **Breadth ≥10x** | % covered at ≥10x depth (confident detection) |
| **Breadth ≥30x** | % covered at ≥30x depth (assembly-ready) |
| **Mean depth** | Average depth across all genome positions |
| **Mapped reads** | Total reads mapping to the organism (all contigs aggregated) |

### Why minimap2, not KMA

KMA (used in step 11 for AMR screening against ResFinder gene sequences) was initially tested for phenotypic verification against whole genomes. The results were unreliable — KMA's competitive mapping and scoring model is optimized for short gene-length references, not complete chromosomes. minimap2 with `map-ont` preset provides real genome-wide coverage metrics (breadth, depth) that are biologically interpretable.

---

## 5. Reports

### EpiTax Report (step 14)
- **Format:** Self-contained HTML (~5 MB, Plotly.js embedded, works offline)
- **Sections:** 9 (survival, contamination, taxonomy, diversity, viability, AMR, phenotypic, rarefaction, conclusion)
- **Features:** Sortable tables, sample filter checkboxes, Whittaker rank-abundance curves, AMR stacked bars by antibiotic class, MAG depth candidates with 30x threshold line
- **CLI flag:** `--results-dir` enables running outside Nextflow for testing

### EpiTaxMAG Report (step 32)
- **HTML:** 16 sections with sample filter, search box, screening charts
- **Excel:** 16 sheets with all raw data (openpyxl)
- **PDF:** Formal report with institutional branding (ReportLab)

---

## 6. Databases

Total size: ~170 GB (TAX only) or ~360 GB (TAX + MAG).

| Database | Size | Used by | Auto-download |
|----------|------|---------|---------------|
| Kraken2 PlusPF (2025-10-15) | ~70 GB | KRAKEN2, BRACKEN | Yes |
| Kaiju RefSeq | ~30 GB | KAIJU | Yes |
| Sylph GTDB r226 + Fungi + Viral | ~15 GB | SYLPH | Yes |
| ResFinder (DTU) | ~3 MB | KMA_READS | Yes |
| FastQ Screen panel (Zenodo) | ~2 GB | FASTQSCREEN | Yes |
| Phenotypic type-strains | ~50 MB | MINIMAP2_PHENOTYPIC | Yes (from NCBI) |
| CheckM2 | ~3 GB | CHECKM2 | Yes |
| GTDB-Tk r232 | ~100 GB | GTDBTK | Yes |
| Sourmash GTDB RS226 | ~3.7 GB | SOURMASH_CLASSIFY | Yes |
| AMRFinderPlus | ~1 GB | AMRFINDERPLUS_MAGS | Yes |
| Bakta full | ~84 GB | BAKTA | Yes |
| geNomad | ~3 GB | GENOMAD | Yes |

Default location: `$projectDir/databases/`. Override with `--db_root /path/to/dbs`.

Automatic download on first run (via `workflow/setup_db.nf`). Disable with `--skip_db_setup`.

---

## 7. Hardware Profiles

| Profile | CPUs | RAM | Kraken2 mode | Use case |
|---------|------|-----|--------------|----------|
| `gru` | 48 | 240 GB | Full RAM | High-performance server |
| `censalud` | 24 | 58 GB | memory-mapping (NVMe) | CenSalud i7-14700 |
| `pc` | 20 | 52 GB | memory-mapping | Desktop workstation |

Select with `-profile gru`. Combine with Singularity: `-profile gru,singularity`.

GTDB-Tk requires >50 GB RAM. On machines with less, use `--skip_gtdbtk` (Sourmash alternative, ~5 GB RAM).

---

## 8. Singularity & Conda

### Conda (default)
Each tool has its own conda environment YAML in `envs/`. Environments are created and cached automatically in `work/conda/` on first run.

### Singularity
Activated with `-profile singularity`. All bioinformatics tools have containers from Galaxy Depot. Python reporting processes (CHECK_SQCORE, FINAL_REPORT) always use conda (no suitable biocontainer for plotly+reportlab+pandas).

A custom image (`epitaxmag-tools`) covers multi-tool processes: FastQ Screen + minimap2 + samtools + KMA + Python reporting. Build with `bash containers/build.sh docker`.

Bind paths and container path configured in `conf/local.config`.

---

## 9. Key Design Decisions

### No Racon polishing
Medaka 2.x models for R10.4.1 SUP v5 are trained on raw Flye output. Adding Racon rounds can degrade quality. The polishing chain is simply: Flye → Medaka.

### Medaka threading
`medaka inference` is single-threaded (neural network). We allocate 2 CPUs per task to maximize parallelism across samples rather than wasting threads.

### Three-binner consensus
MetaBAT2 (composition+coverage) + MaxBin2 (EM) + SemiBin2 (deep learning, long-read mode) → DAS Tool consensus. Each binner has different strengths; the consensus improves bin quality.

### SemiBin2 long-read mode
`--sequencing-type long_read` is critical. Without it, SemiBin2 uses short-read defaults that produce poor bins from Nanopore data.

### GUNC as post-filter
GUNC (chimerism detection) acts as a post-filter: bins flagged as chimeric are discarded even if they pass MIMAG quality thresholds.

### KMA for AMR, minimap2 for phenotypic
KMA with `-ont` mode is optimal for detecting short gene sequences (AMR, ResFinder). For whole-genome verification (phenotypic), minimap2 provides meaningful breadth/depth metrics across complete chromosomes.

### Conditional steps
- Host removal: only if `--host_genome` is set
- Phenotypic verification: on by default, disable with `--phenotypic_targets false`
- GTDB-Tk vs Sourmash: controlled by `--skip_gtdbtk`
- MAG phase: only with `--run_assembly`

---

## 10. Process Reference Table

Complete mapping of Nextflow process names → output directories → tools.

### Phase TAX

| Step | Process Name | Directory | Tools |
|------|-------------|-----------|-------|
| 01 | `EPITAXMAG_TAX:NANOPLOT_RAW` | `01_qc_raw/` | NanoPlot 1.46.2 |
| 01 | `EPITAXMAG_TAX:FASTQC_RAW` | `01_qc_raw/` | FastQC 0.12.1 |
| 01 | `EPITAXMAG_TAX:CHECK_SQCORE_RAW` | `01_qc_raw/` | Python |
| 02 | `EPITAXMAG_TAX:FASTQSCREEN` | `02_fastqscreen_raw/` | FastQ Screen 0.16.0 |
| 03 | `EPITAXMAG_TAX:PORECHOP` | `03_trim_porechop/` | Porechop ABI 0.5.1 |
| 04 | `EPITAXMAG_TAX:CHOPPER` | `04_filter_chopper/` | Chopper 0.12.0 |
| 05 | `EPITAXMAG_TAX:HOST_REMOVAL` | `05_host_removal/` | minimap2 + samtools |
| 06 | `EPITAXMAG_TAX:NANOPLOT_FILTERED` | `06_qc_filtered/` | NanoPlot 1.46.2 |
| 06 | `EPITAXMAG_TAX:FASTQC_FILTERED` | `06_qc_filtered/` | FastQC 0.12.1 |
| 06 | `EPITAXMAG_TAX:CHECK_SQCORE_FILTERED` | `06_qc_filtered/` | Python |
| 07 | `EPITAXMAG_TAX:KRAKEN2` | `07_tax_kraken2/` | Kraken2 2.17.1 |
| 08 | `EPITAXMAG_TAX:BRACKEN` | `08_tax_bracken/` | Bracken 3.1 |
| 09 | `EPITAXMAG_TAX:KAIJU` | `09_tax_kaiju/` | Kaiju 1.10.1 |
| 10 | `EPITAXMAG_TAX:SYLPH` | `10_tax_sylph/` | Sylph 0.9.0 |
| 11 | `EPITAXMAG_TAX:KMA_READS` | `11_amr_kma/` | KMA 1.6.8 |
| 12 | `EPITAXMAG_TAX:BUILD_PHENOTYPIC_DB` | `12_phenotypic_minimap2/db/` | Python + NCBI API |
| 12 | `EPITAXMAG_TAX:MINIMAP2_PHENOTYPIC` | `12_phenotypic_minimap2/` | minimap2 + samtools |
| 13 | `EPITAXMAG_TAX:MULTIQC` | `13_multiqc/` | MultiQC 1.33 |
| 14 | `EPITAXMAG_TAX:FINAL_REPORT` | `14_epitax_report/` | Python (Plotly) |

### Phase MAG

| Step | Process Name | Directory | Tools |
|------|-------------|-----------|-------|
| 15 | `EPITAXMAG_MAG:METAFLYE` | `15_assembly_flye/` | Flye 2.9.6 |
| 16 | `EPITAXMAG_MAG:MEDAKA` | `16_polish_medaka/` | Medaka 2.2.1 |
| 17 | `EPITAXMAG_MAG:QUAST` | `17_qc_assembly/` | QUAST 5.3.0 |
| 18 | `EPITAXMAG_MAG:MAP_READS` | `18_coverage/` | minimap2 + samtools |
| 19 | `EPITAXMAG_MAG:METABAT2` | `19_binning_metabat2/` | MetaBAT2 2.18 |
| 20 | `EPITAXMAG_MAG:MAXBIN2` | `20_binning_maxbin2/` | MaxBin2 2.2.7 |
| 21 | `EPITAXMAG_MAG:SEMIBIN2` | `21_binning_semibin2/` | SemiBin2 2.2.1 |
| 22 | `EPITAXMAG_MAG:DAS_TOOL` | `22_binning_dastool/` | DAS Tool 1.1.7 |
| 23 | `EPITAXMAG_MAG:CHECKM2` | `23_binqc_checkm2/` | CheckM2 1.1.0 |
| 24 | `EPITAXMAG_MAG:GTDBTK` | `24_taxonomy_gtdbtk/` | GTDB-Tk 2.7.0 |
| 24 | `EPITAXMAG_MAG:SOURMASH_CLASSIFY` | `24_taxonomy_sourmash/` | Sourmash 4.9.4 |
| 25 | `EPITAXMAG_MAG:AMRFINDERPLUS_MAGS` | `25_amr_mags/` | AMRFinderPlus 4.2.7 |
| 26 | `EPITAXMAG_MAG:GENOMAD` | `26_genomad/` | geNomad 1.12.0 |
| 27 | `EPITAXMAG_MAG:AMR_PATHOGEN_INTEGRATION` | `27_amr_pathogen_integration/` | Python |
| 28 | `EPITAXMAG_MAG:ABRICATE` | `28_abricate/` | ABRicate 1.4.0 |
| 29 | `EPITAXMAG_MAG:MOBSUITE` | `29_mobsuite/` | MOB-suite 3.1.9 |
| 30 | `EPITAXMAG_MAG:INTEGRONFINDER` | `30_integronfinder/` | IntegronFinder 2.0 |
| 31 | `EPITAXMAG_MAG:BAKTA` | `31_bakta/` | Bakta 1.12.0 |
| 32 | `EPITAXMAG_MAG:MAG_REPORT` | `32_reports/` | Python (Plotly+Excel+PDF) |
