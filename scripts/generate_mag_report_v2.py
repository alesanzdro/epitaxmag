#!/usr/bin/env python3
"""
generate_mag_report_v2.py
Generates a comprehensive, self-contained, interactive HTML report for the
EpiTaxMAG pipeline.  Single-file output with embedded Plotly.js, interactive
tables (sortable + filterable), and professional EpiTax branding.

Usage:
    python3 scripts/generate_mag_report_v2.py \
        --results-dir results/260226_EPIM232 \
        --run-name 260226_EPIM232 \
        --input-dir /home/asanzc/jobs/nanometa/data/260226_EPIM232 \
        --branding branding/ \
        --output results/260226_EPIM232/32_reports/260226_EPIM232_epitaxmag_report_v2.html
"""

import argparse, os, sys, glob, base64, warnings, re
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly

# Import screening charts module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_screening_charts import (load_fastqscreen_raw, build_fastqscreen_100pct,
                                        load_sylph_profile, build_sylph_dual_heatmap)

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════
# PALETTE & CONSTANTS
# ═══════════════════════════════════════════════════════════════

BLUE = '#2C5F8A'; LBLUE = '#5B9ABF'; ORANGE = '#E05C2A'
GREEN = '#28A745'; RED = '#E74C3C'; GREY = '#6C757D'
PALETTE = [BLUE, ORANGE, GREEN, '#9B59B6', '#F39C12',
           '#1ABC9C', RED, '#3498DB', '#95A5A6', '#D35400',
           '#27AE60', '#8E44AD', '#C0392B', '#16A085', '#2980B9']

RISK_COLORS = {'CRITICAL': RED, 'HIGH': ORANGE, 'MEDIUM': '#F39C12', 'LOW': GREEN}
RISK_BADGE = {
    'CRITICAL': f'<span class="badge" style="background:{RED}">CRITICAL</span>',
    'HIGH':     f'<span class="badge" style="background:{ORANGE}">HIGH</span>',
    'MEDIUM':   f'<span class="badge" style="background:#F39C12">MEDIUM</span>',
    'LOW':      f'<span class="badge" style="background:{GREEN}">LOW</span>',
}
QC_BADGE = {
    'HQ': f'<span class="badge" style="background:{GREEN}">HQ</span>',
    'MQ': f'<span class="badge" style="background:{ORANGE}">MQ</span>',
    'LQ': f'<span class="badge" style="background:{GREY}">LQ</span>',
}
CONCORDANCE_BADGE = {
    'BOTH':    f'<span class="badge" style="background:{GREEN}">Reads+Contigs</span>',
    'READS':   f'<span class="badge" style="background:#F39C12">Reads only</span>',
    'CONTIGS': f'<span class="badge" style="background:{LBLUE}">Contigs only</span>',
}

# ═══════════════════════════════════════════════════════════════
# PROVENANCE — single source of truth for "which tool produced what"
# Each entry knows: pinned version, role, reference URL, layer
# (reads / contigs / mags / cross-cutting), and which Nextflow params
# materially change its result. The HTML provenance panels render
# directly from this registry, so adding a new tool == one edit here.
# ═══════════════════════════════════════════════════════════════

TOOL_REGISTRY = {
    # key                : (label, version, role, layer, url, params, process_name)
    # `process_name` is the Nextflow process name (suffix only, without the
    # workflow prefix) used to join with trace.txt. None when the tool runs
    # as part of a multi-tool process and has no dedicated entry in trace.
    'nanoplot'           : ('NanoPlot',        '1.46.2',     'Read QC (length, quality)',              'reads',      'https://github.com/wdecoster/NanoPlot',                                          [],                                       'NANOPLOT_RAW|NANOPLOT_FILTERED'),
    'fastqc'             : ('FastQC',          '0.12.1',     'Read QC (per-base quality, GC, dups)',   'reads',      'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/',                     [],                                       'FASTQC_RAW|FASTQC_FILTERED'),
    'fastqscreen'        : ('FastQ Screen',    '0.16.0',     'Cross-contamination screening',          'reads',      'https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/',               ['fastqscreen_conf'],                     'FASTQSCREEN'),
    'porechop'           : ('Porechop ABI',    '0.5.1',      'Adapter trimming (ONT)',                  'reads',      'https://github.com/bonsai-team/Porechop_ABI',                                    [],                                       'PORECHOP'),
    'chopper'            : ('Chopper',         '0.12.0',     'Quality and length filtering',            'reads',      'https://github.com/wdecoster/chopper',                                            ['min_quality', 'min_length'],            'CHOPPER'),
    'minimap2_host'      : ('minimap2',        '2.30',       'Host removal (read mapping)',             'reads',      'https://github.com/lh3/minimap2',                                                ['host_genome'],                          'HOST_REMOVAL'),
    'samtools'           : ('samtools',        '1.21',       'BAM/SAM manipulation',                    'cross',      'https://www.htslib.org/',                                                        [],                                       None),
    'kraken2'            : ('Kraken2',         '2.17.1',     'Read taxonomic classification (k-mer)',   'reads',      'https://ccb.jhu.edu/software/kraken2/',                                          ['kraken2_db', 'kraken2_confidence'],      'KRAKEN2'),
    'bracken'            : ('Bracken',         '3.1',        'Abundance re-estimation on Kraken2',      'reads',      'https://ccb.jhu.edu/software/bracken/',                                          ['kraken2_db'],                           'BRACKEN'),
    'kaiju'              : ('Kaiju',           '1.10.1',     'Protein-level read taxonomy',             'reads',      'https://bioinformatics-centre.github.io/kaiju/',                                 ['kaiju_nodes', 'kaiju_names', 'kaiju_fmi'], 'KAIJU'),
    'sylph'              : ('Sylph',           '0.9.0',      'k-mer profiling and ANI estimation',      'reads',      'https://github.com/bluenote-1577/sylph',                                         ['sylph_db', 'sylph_db_fungi', 'sylph_db_viral'], 'SYLPH'),
    'kma'                : ('KMA',             '1.6.8',      'AMR detection in reads (ResFinder)',      'reads',      'https://bitbucket.org/genomicepidemiology/kma/',                                 ['resfinder_fsa'],                        'KMA_READS'),
    'multiqc'            : ('MultiQC',         '1.33',       'QC aggregation',                          'reads',      'https://multiqc.info/',                                                          [],                                       'MULTIQC'),
    'minimap2_phenotypic': ('minimap2',        '2.30',       'Phenotypic verification (whole-genome)',  'reads',      'https://github.com/lh3/minimap2',                                                ['phenotypic_targets'],                   'MINIMAP2_PHENOTYPIC|BUILD_PHENOTYPIC_DB'),
    'metaflye'           : ('MetaFlye',        '2.9.6',      'Metagenomic assembly',                    'contigs',    'https://github.com/fenderglass/Flye',                                            ['flye_mode', 'genome_size_mb', 'min_coverage_mag'], 'METAFLYE'),
    'medaka'             : ('Medaka',          '2.2.1',      'Consensus polishing (ONT)',               'contigs',    'https://github.com/nanoporetech/medaka',                                          ['medaka_model'],                         'MEDAKA'),
    'quast'              : ('QUAST',           '5.3.0',      'Assembly metrics',                        'contigs',    'https://quast.sourceforge.net/',                                                  [],                                       'QUAST'),
    'metabat2'           : ('MetaBAT2',        '2.18',       'Coverage-based binning',                  'mags',       'https://bitbucket.org/berkeleylab/metabat/',                                      [],                                       'METABAT2'),
    'maxbin2'            : ('MaxBin2',         '2.2.7',      'EM-based binning',                        'mags',       'https://sourceforge.net/projects/maxbin2/',                                       [],                                       'MAXBIN2'),
    'semibin2'           : ('SemiBin2',        '2.2.1',      'Deep-learning binning (long reads)',      'mags',       'https://github.com/BigDataBiology/SemiBin',                                       [],                                       'SEMIBIN2'),
    'dastool'            : ('DAS Tool',        '1.1.7',      'Bin refinement (consensus)',              'mags',       'https://github.com/cmks/DAS_Tool',                                                [],                                       'DAS_TOOL'),
    'checkm2'            : ('CheckM2',         '1.1.0',      'MAG completeness and contamination',      'mags',       'https://github.com/chklovski/CheckM2',                                            ['checkm2_db', 'mag_completeness_min', 'mag_contamination_max'], 'CHECKM2'),
    'skani'              : ('Skani',           '0.3.1',      'ANI-based MAG taxonomy (GTDB)',           'mags',       'https://github.com/bluenote-1577/skani',                                          ['taxonomy_tool', 'skani_db'],            'SKANI_CLASSIFY'),
    'sourmash'           : ('Sourmash',        '4.9.5',      'k-mer sketches MAG taxonomy (GTDB)',      'mags',       'https://github.com/sourmash-bio/sourmash',                                        ['taxonomy_tool'],                        'SOURMASH_CLASSIFY'),
    'gtdbtk'             : ('GTDB-Tk',         '2.7.0',      'Phylogenetic MAG taxonomy (GTDB)',        'mags',       'https://ecogenomics.github.io/GTDBTk/',                                           ['taxonomy_tool', 'gtdbtk_data'],         'GTDBTK'),
    'amrfinderplus'      : ('AMRFinderPlus',   '4.2.7',      'AMR genes on contigs (NCBI)',             'mags',       'https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/',      ['amrfinder_db'],                         'AMRFINDERPLUS_MAGS|AMR_ON_PLASMIDS'),
    'genomad'            : ('geNomad',         '1.12.0',     'Plasmid and virus prediction',            'contigs',    'https://github.com/apcamargo/genomad',                                            ['genomad_db'],                           'GENOMAD'),
    'abricate'           : ('ABRicate',        '1.0.1',      'Virulence (VFDB) and AMR (CARD)',         'mags',       'https://github.com/tseemann/abricate',                                            [],                                       'ABRICATE'),
    'mobsuite'           : ('MOB-suite',       '3.1.9',      'Plasmid typing (replicon, mobility)',     'mags',       'https://github.com/phac-nml/mob-suite',                                          [],                                       'MOBSUITE'),
    'integronfinder'     : ('IntegronFinder',  '2.0.5',      'Integron and gene cassette detection',    'mags',       'https://github.com/gem-pasteur/Integron_Finder',                                  [],                                       'INTEGRONFINDER'),
    'bakta'              : ('Bakta',           '1.12.0',     'MAG functional annotation',               'mags',       'https://github.com/oschwengers/bakta',                                            ['bakta_db', 'bakta_db_type'],            'BAKTA'),
}

LAYER_LABELS = {
    'reads'  : 'Reads layer',
    'contigs': 'Contigs / assembly layer',
    'mags'   : 'MAGs layer',
    'cross'  : 'Cross-cutting',
}

# Per-section provenance: which tools produced the data shown.
# Section keys must match the section's anchor IDs in the HTML.
SECTION_TOOLS = {
    'qc-survival'  : ['nanoplot', 'fastqc', 'porechop', 'chopper', 'minimap2_host'],
    'screening'    : ['fastqscreen', 'sylph'],
    'fastqscreen'  : ['fastqscreen'],
    'sylph'        : ['sylph'],
    'reads-tax'    : ['kraken2', 'bracken', 'kaiju'],
    'phenotypic'   : ['minimap2_phenotypic', 'samtools'],
    'amr-reads'    : ['kma'],
    'assembly'     : ['metaflye', 'medaka', 'quast'],
    'binning'      : ['metabat2', 'maxbin2', 'semibin2', 'dastool'],
    'mags'         : ['checkm2', 'bakta'],
    'mag-taxonomy' : ['skani', 'sourmash', 'gtdbtk'],
    'amr-mags'     : ['amrfinderplus'],
    'amr-heatmap'  : ['amrfinderplus'],
    'concordance'  : ['kma', 'amrfinderplus'],
    'vfdb'         : ['abricate'],
    'plasmids'     : ['genomad', 'amrfinderplus'],
    'mobsuite'     : ['mobsuite'],
    'integrons'    : ['integronfinder'],
    'risk'         : ['amrfinderplus', 'genomad', 'checkm2'],
    'multiqc'      : ['multiqc'],
}


def _tool_entry_for_software_table(tool_key, params_payload=None):
    """Return (label, version_str, role, url) tuple, with parameter
    values appended to the version string when params_payload is provided."""
    if tool_key not in TOOL_REGISTRY:
        return None
    label, version, role, _layer, url, _, _ = TOOL_REGISTRY[tool_key]
    return label, version, role, url


def _parse_duration(s):
    """Convert Nextflow duration strings ('1m 32s', '20.1s', '2h 15m')
    to seconds. Returns 0.0 on failure."""
    if s is None:
        return 0.0
    s = str(s).strip()
    if not s or s == '-':
        return 0.0
    try:
        return float(s)
    except ValueError:
        pass
    total = 0.0
    # `ms` must precede `m` in the alternation so "649ms" is not parsed
    # as "649 m". Same trick for the other multi-letter units.
    parts = re.findall(r'(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)', s)
    for value, unit in parts:
        v = float(value)
        if unit == 'd':
            total += v * 86400
        elif unit == 'h':
            total += v * 3600
        elif unit == 'm':
            total += v * 60
        elif unit == 's':
            total += v
        elif unit == 'ms':
            total += v / 1000.0
    return total


def _parse_size_bytes(s):
    """Convert Nextflow memory strings ('482.4 MB', '2.1 GB') to bytes."""
    if s is None:
        return 0
    s = str(s).strip()
    if not s or s == '-':
        return 0
    try:
        return int(float(s))
    except ValueError:
        pass
    m = re.match(r'(\d+(?:\.\d+)?)\s*([KMGTP]?B)', s, re.IGNORECASE)
    if not m:
        return 0
    value = float(m.group(1))
    unit = m.group(2).upper()
    factors = {'B': 1, 'KB': 1024, 'MB': 1024**2,
               'GB': 1024**3, 'TB': 1024**4, 'PB': 1024**5}
    return int(value * factors.get(unit, 1))


def load_trace(results_dir):
    """Parse Nextflow trace.txt into a DataFrame with normalised seconds
    and bytes columns. Returns empty df when file is missing."""
    candidates = [
        f'{results_dir}/00_reports/trace.txt',
        f'{results_dir}/trace.txt',
    ]
    path = next((p for p in candidates if os.path.isfile(p)), None)
    if path is None:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep='\t')
    except Exception as e:
        print(f'[WARN] Could not parse trace.txt: {e}')
        return pd.DataFrame()
    if df.empty:
        return df

    df['Process_short'] = df['process'].astype(str).apply(
        lambda v: v.split(':', 1)[1] if ':' in v else v)
    df['Process_layer'] = df['process'].astype(str).apply(
        lambda v: (v.split(':', 1)[0] if ':' in v else 'OTHER').replace(
            'EPITAXMAG_', '') or 'OTHER')

    df['Duration_s'] = df.get('duration', '').apply(_parse_duration)
    df['Realtime_s'] = df.get('realtime', '').apply(_parse_duration)
    df['Peak_rss_b'] = df.get('peak_rss', '').apply(_parse_size_bytes)
    df['Rss_b']      = df.get('rss', '').apply(_parse_size_bytes)
    return df


# Layer presentation order for the Pipeline-execution section
TRACE_LAYERS = ['SETUP', 'TAX', 'MAG']

# Anomaly thresholds — surface real outliers on the validation runs
# (Kaiju RAM, GTDB-Tk when enabled, slow Bakta bins) without flagging
# the ordinary heavy-but-expected processes.
TRACE_ANOMALY_RAM_GB = 40
TRACE_ANOMALY_DURATION_S = 30 * 60


def trace_summary(trace_df):
    if trace_df.empty:
        return {}
    n_total = len(trace_df)
    n_cached = int((trace_df['status'] == 'CACHED').sum())
    n_completed = int((trace_df['status'] == 'COMPLETED').sum())
    n_failed = int((trace_df['status'].astype(str)
                    .str.upper().isin({'FAILED', 'ERROR'})).sum())
    pct_cached = (100 * n_cached / n_total) if n_total else 0
    cpu_hours = float(trace_df['Realtime_s'].sum()) / 3600.0
    return {
        'n_total': n_total,
        'n_cached': n_cached,
        'n_completed': n_completed,
        'n_failed': n_failed,
        'pct_cached': pct_cached,
        'cpu_hours': cpu_hours,
    }


def trace_top_n(trace_df, by, n=5, status_filter=None):
    if trace_df.empty:
        return pd.DataFrame()
    df = trace_df.copy()
    if status_filter is not None:
        df = df[df['status'].isin(status_filter)]
    if df.empty:
        return df
    return df.nlargest(n, by)


def trace_anomalies(trace_df, p_quantile=0.95):
    """Combine absolute thresholds (always-flag) with run-relative
    detection (top-N% within this run). The relative quantile is computed
    only over actually-executed tasks (status != CACHED would mean we
    miss the historical CACHED rows that hold the original measurements
    — so we keep both, but compute the quantile over the union).

    Returns df with extra column `Anomaly_flags` (list of strings:
    RAM, SLOW, RAM-P95, SLOW-P95)."""
    if trace_df.empty:
        return pd.DataFrame()

    df = trace_df.copy()
    df['Anomaly_flags'] = [[] for _ in range(len(df))]

    ram_abs = TRACE_ANOMALY_RAM_GB * (1024 ** 3)
    df.loc[df['Peak_rss_b'] > ram_abs, 'Anomaly_flags'] = df.loc[
        df['Peak_rss_b'] > ram_abs, 'Anomaly_flags'].apply(lambda l: l + ['RAM'])
    df.loc[df['Realtime_s'] > TRACE_ANOMALY_DURATION_S, 'Anomaly_flags'] = df.loc[
        df['Realtime_s'] > TRACE_ANOMALY_DURATION_S, 'Anomaly_flags'].apply(
            lambda l: l + ['SLOW'])

    # Run-relative: only meaningful with enough rows. Skip near-zero
    # values so a one-shot 30-second helper task does not dominate the
    # P95 of a run dominated by CACHED entries.
    nonzero_rt = df[df['Realtime_s'] > 1.0]['Realtime_s']
    if len(nonzero_rt) >= 10:
        q_rt = nonzero_rt.quantile(p_quantile)
        df.loc[df['Realtime_s'] > q_rt, 'Anomaly_flags'] = df.loc[
            df['Realtime_s'] > q_rt, 'Anomaly_flags'].apply(
                lambda l: l + [f'SLOW-P{int(p_quantile * 100)}'])
    nonzero_rss = df[df['Peak_rss_b'] > 0]['Peak_rss_b']
    if len(nonzero_rss) >= 10:
        q_rss = nonzero_rss.quantile(p_quantile)
        df.loc[df['Peak_rss_b'] > q_rss, 'Anomaly_flags'] = df.loc[
            df['Peak_rss_b'] > q_rss, 'Anomaly_flags'].apply(
                lambda l: l + [f'RAM-P{int(p_quantile * 100)}'])

    return df[df['Anomaly_flags'].apply(len) > 0].copy()


def trace_for_tool(trace_df, tool_key):
    """Rows from trace_df belonging to a tool, via the pipe-separated
    process name in TOOL_REGISTRY."""
    if trace_df.empty or tool_key not in TOOL_REGISTRY:
        return pd.DataFrame()
    proc_name = TOOL_REGISTRY[tool_key][6]
    if not proc_name:
        return pd.DataFrame()
    names = set(proc_name.split('|'))
    return trace_df[trace_df['Process_short'].isin(names)]


def _fmt_seconds(seconds):
    s = float(seconds or 0)
    if s < 60:
        return f'{s:.1f}s'
    if s < 3600:
        return f'{int(s // 60)}m {int(s % 60)}s'
    h = int(s // 3600)
    m = int((s % 3600) // 60)
    return f'{h}h {m}m'


def _fmt_bytes(b):
    b = float(b or 0)
    for unit, label in [(1024**4, 'TB'), (1024**3, 'GB'),
                        (1024**2, 'MB'), (1024, 'KB')]:
        if b >= unit:
            return f'{b / unit:.1f} {label}'
    return f'{int(b)} B'


# ─── FastQ Screen detail (interpreted) ─────────────────────────
# We classify each genome of the FastQ Screen panel into a category so
# the report can separate "this is what we expected to see" (surveillance
# target organisms) from "this is contamination we need to flag" (host,
# technical controls, vectors). The classification IS NOT in the
# *_screen.txt — it is an inference made by this report from the
# `params.target_organisms` list and a small known-panel lookup.

# Common short names from the curated panel → genus name (so we can match
# against `target_organisms`). Add new entries as the panel grows.
_PANEL_GENUS_LOOKUP = {
    'A_baumannii'        : 'Acinetobacter',
    'A_bouvetii'         : 'Acinetobacter',
    'Aeromonas'          : 'Aeromonas',
    'Aliarcobacter'      : 'Aliarcobacter',
    'Arcobacter_butzleri': 'Arcobacter',
    'B_graminisolvens'   : 'Bacteroides',
    'B_ovatus'           : 'Bacteroides',
    'B_pseudomallei'     : 'Burkholderia',
    'Campylobacter'      : 'Campylobacter',
    'Dorea'              : 'Dorea',
    'Ecoli'              : 'Escherichia',
    'E_faecalis'         : 'Enterococcus',
    'E_faecium'          : 'Enterococcus',
    'H_sapiens'          : 'Homo',
    'Klebsiella'         : 'Klebsiella',
    'Legionella'         : 'Legionella',
    'Leptospira'         : 'Leptospira',
    'M_avium'            : 'Mycobacterium',
    'N_fowleri'          : 'Naegleria',
    'P_aeruginosa'       : 'Pseudomonas',
    'Phocaeicola_dorei'  : 'Phocaeicola',
    'Salmonella'         : 'Salmonella',
    'Shigella'           : 'Shigella',
    'Vibrio_spp'         : 'Vibrio',
}

_HOST_PANEL = {'H_sapiens', 'Mus_musculus', 'Bos_taurus', 'Sus_scrofa', 'Gallus'}
_CONTROL_PANEL = {'Controls', 'PhiX', 'Avian_Marker', 'Vectors'}
_PARASITE_PANEL = {'N_fowleri', 'Parasites', 'Cryptosporidium', 'Giardia'}
_VIRUS_PANEL = {'Adenovirus', 'Virus_Enteric'}


def classify_panel_genome(name, target_genera):
    """Return (category, genus_or_label). Categories:
    'target' (surveillance), 'host', 'control', 'parasite', 'virus',
    'other_reference'."""
    if name in _CONTROL_PANEL:
        return 'control', name
    if name in _HOST_PANEL:
        return 'host', name
    if name in _PARASITE_PANEL:
        return 'parasite', name
    if name in _VIRUS_PANEL:
        return 'virus', name
    inferred_genus = _PANEL_GENUS_LOOKUP.get(name) or name.replace('_', ' ').split()[0]
    if any(g.lower() == inferred_genus.lower() or
           g.lower() in name.lower()
           for g in target_genera):
        return 'target', inferred_genus
    return 'other_reference', inferred_genus


def load_bracken_genus(results_dir):
    """Load per-sample genus-level Bracken counts. Returns long DataFrame
    with Sample · Genus · Reads · Pct."""
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/08_tax_bracken/*.bracken.G.txt')):
        sample = os.path.basename(f).replace('.bracken.G.txt', '')
        try:
            df = pd.read_csv(f, sep='\t')
        except Exception:
            continue
        if df.empty or 'name' not in df.columns:
            continue
        for _, r in df.iterrows():
            rows.append({
                'Sample': sample,
                'Genus' : str(r['name']).strip(),
                'Reads' : int(r.get('new_est_reads', 0)),
                'Pct'   : float(r.get('fraction_total_reads', 0)) * 100,
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


_SYLPH_GENUS_RE = re.compile(
    r'^(?:.*?MAG:\s*)?([A-Z][a-zA-Z0-9_-]+)\b')


def _extract_sylph_genus(contig_name):
    """Best-effort genus extraction from the Sylph Contig_name field.
    Tries 'MAG: Acinetobacter sp. ...' patterns and falls back to the
    first capitalised token. Returns empty string when nothing fits."""
    if not contig_name:
        return ''
    s = str(contig_name)
    if s.startswith('NZ_') or s.startswith('NC_') or s.startswith('GCF_'):
        # The accession itself is at the start; the genus is later.
        s = s.split(' ', 1)[-1] if ' ' in s else s
    m = _SYLPH_GENUS_RE.match(s)
    if m:
        candidate = m.group(1)
        # Reject obvious non-genera (single-letter prefixes, NA, etc.)
        if candidate.upper() in {'NA', 'MAG', 'TPA', 'UNCULTURED'} or len(candidate) < 4:
            return ''
        return candidate
    return ''


def load_sylph_genus(results_dir):
    """Load per-sample genus-level Sylph abundances by aggregating the
    full profile rows. Genus is parsed from the Contig_name field
    (Sylph reports the best-hit reference genome, not a NCBI taxid)."""
    rows = []
    candidate_paths = [
        f'{results_dir}/10_tax_sylph/sylph_profile_all.tsv',
    ]
    candidate_paths += sorted(
        glob.glob(f'{results_dir}/10_tax_sylph/*.sylph.tsv'))
    seen = set()
    for path in candidate_paths:
        if not os.path.isfile(path) or path in seen:
            continue
        seen.add(path)
        try:
            df = pd.read_csv(path, sep='\t')
        except Exception:
            continue
        if df.empty or 'Sample_file' not in df.columns:
            continue
        df['Genus'] = df.get('Contig_name', '').apply(_extract_sylph_genus)
        for (sample_file, genus), grp in df.groupby(['Sample_file', 'Genus']):
            if not genus:
                continue
            sample = os.path.basename(str(sample_file))
            for ext in ('.filtered.fastq.gz', '.fastq.gz', '.fq.gz'):
                if sample.endswith(ext):
                    sample = sample[:-len(ext)]
                    break
            rows.append({
                'Sample': sample,
                'Genus' : genus,
                'Pct'   : float(grp['Taxonomic_abundance'].sum()),
                'Hits'  : int(len(grp)),
                'Mean_ANI': float(pd.to_numeric(
                    grp.get('Adjusted_ANI', pd.Series([])),
                    errors='coerce').mean() or 0),
            })
        # We use the first source that yielded data; the rest would duplicate.
        if rows:
            break
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def compare_bracken_sylph(bracken_df, sylph_df, top_n=15):
    """Per-sample concordance between Bracken and Sylph at genus level.
    Returns a long DataFrame: Sample · Genus · Bracken_pct · Sylph_pct ·
    Detected_by ∈ {Both, Bracken only, Sylph only}."""
    if bracken_df.empty and sylph_df.empty:
        return pd.DataFrame()
    samples = sorted(set(
        list(bracken_df.get('Sample', pd.Series(dtype=str)).unique()) +
        list(sylph_df.get('Sample', pd.Series(dtype=str)).unique())))
    rows = []
    for sample in samples:
        b = bracken_df[bracken_df['Sample'] == sample]
        s = sylph_df[sylph_df['Sample'] == sample]
        # Restrict to top-N most abundant on each side, union those keys.
        b_top = set(b.nlargest(top_n, 'Pct')['Genus']) if not b.empty else set()
        s_top = set(s.nlargest(top_n, 'Pct')['Genus']) if not s.empty else set()
        union = sorted(b_top | s_top)
        b_lookup = b.set_index('Genus')['Pct'].to_dict() if not b.empty else {}
        s_lookup = s.set_index('Genus')['Pct'].to_dict() if not s.empty else {}
        for genus in union:
            in_b = genus in b_top
            in_s = genus in s_top
            detected = ('Both' if in_b and in_s
                        else ('Bracken only' if in_b else 'Sylph only'))
            rows.append({
                'Sample'      : sample,
                'Genus'       : genus,
                'Bracken_pct' : b_lookup.get(genus, 0.0),
                'Sylph_pct'   : s_lookup.get(genus, 0.0),
                'Detected_by' : detected,
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_fastqscreen_detail(results_dir):
    """Parse every *_screen.txt under 02_fastqscreen_raw/ into a tidy
    long DataFrame with one row per (sample, genome) and the percentages
    already split into 5 mutually-exclusive categories straight from
    the file. Header columns from FastQ Screen 0.16.0.

    Returns a tuple (detail_df, no_hit_pct_per_sample) where the second
    is a dict {sample: %Hit_no_genomes from the file footer}.

    Columns produced (all percentages, 0–100):
      Sample · Genome · Reads_processed · Unmapped_pct ·
      Unique_pct · MultiInternal_pct · CrossMapping_pct · RepeatLike_pct ·
      Total_hit_pct (derived: 100 - Unmapped_pct)
    """
    rows = []
    no_hit = {}
    for f in sorted(glob.glob(f'{results_dir}/02_fastqscreen_raw/*_screen.txt')):
        sample = os.path.basename(f).replace('_screen.txt', '')
        try:
            with open(f) as fh:
                lines = [l.rstrip('\n') for l in fh.readlines()]
        except Exception:
            continue
        # Find header line (starts with 'Genome\t')
        header_idx = None
        for i, l in enumerate(lines):
            if l.startswith('Genome\t'):
                header_idx = i
                break
        if header_idx is None:
            continue
        for l in lines[header_idx + 1:]:
            l = l.strip()
            if not l:
                continue
            if l.startswith('%Hit_no_genomes'):
                m = re.search(r'%Hit_no_genomes:\s*([\d.]+)', l)
                if m:
                    try:
                        no_hit[sample] = float(m.group(1))
                    except ValueError:
                        pass
                continue
            if l.startswith('#') or l.startswith('%'):
                continue
            parts = l.split('\t')
            if len(parts) < 12:
                continue
            try:
                rows.append({
                    'Sample'           : sample,
                    'Genome'           : parts[0],
                    'Reads_processed'  : int(parts[1]),
                    'Unmapped_pct'     : float(parts[3]),
                    'Unique_pct'       : float(parts[5]),
                    'MultiInternal_pct': float(parts[7]),
                    'CrossMapping_pct' : float(parts[9]),
                    'RepeatLike_pct'   : float(parts[11]),
                    'Total_hit_pct'    : 100.0 - float(parts[3]),
                })
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    if not df.empty:
        df.attrs['no_hit_pct'] = no_hit
    return df


def _hurlbert_rarefaction(counts, n_subsample):
    """Expected number of distinct taxa in a random subsample of size n,
    given the per-taxon counts (Hurlbert 1971). Closed-form, no Monte
    Carlo needed.

      E[S(n)] = Σ_i [ 1 − C(N − N_i, n) / C(N, n) ]

    where N is the total reads and N_i is the count of taxon i.
    Implemented with log-binomial coefficients to avoid overflow on
    large samples (N ~ 10^6+).
    """
    counts = np.asarray([c for c in counts if c > 0], dtype=float)
    N = float(counts.sum())
    if N <= 0 or n_subsample <= 0:
        return 0.0
    n = float(min(n_subsample, N))
    if n >= N:
        return float(len(counts))

    from math import lgamma
    def lcomb(a, b):
        if b < 0 or b > a:
            return float('-inf')
        return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)

    log_total = lcomb(N, n)
    if not np.isfinite(log_total):
        return float(len(counts))
    contribs = []
    for ni in counts:
        log_miss = lcomb(N - ni, n) - log_total
        # Numerical floor — extremely abundant taxa always present
        miss = float(np.exp(log_miss)) if np.isfinite(log_miss) else 0.0
        contribs.append(1.0 - min(max(miss, 0.0), 1.0))
    return float(sum(contribs))


def compute_rarefaction(bracken_df, n_steps=20):
    """Per-sample rarefaction curve from Bracken counts.

    Returns long DataFrame: Sample · N_reads · Expected_genera.
    Steps are log-spaced so the early portion of the curve (where the
    inflection lives) is well sampled."""
    if bracken_df is None or bracken_df.empty:
        return pd.DataFrame()
    rows = []
    for sample in sorted(bracken_df['Sample'].unique()):
        sub = bracken_df[bracken_df['Sample'] == sample]
        counts = sub['Reads'].astype(int).tolist()
        total = int(sum(counts))
        if total < 100:
            continue
        steps = np.unique(np.round(
            np.logspace(np.log10(100), np.log10(total), n_steps)).astype(int))
        for n in steps:
            rows.append({
                'Sample': sample,
                'N_reads': int(n),
                'Expected_genera': _hurlbert_rarefaction(counts, n),
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def chart_rarefaction(rarefaction_df):
    if rarefaction_df is None or rarefaction_df.empty:
        return None
    fig = go.Figure()
    samples = sorted(rarefaction_df['Sample'].unique())
    for i, sample in enumerate(samples):
        sub = rarefaction_df[rarefaction_df['Sample'] == sample]
        color = PALETTE[i % len(PALETTE)]
        fig.add_trace(go.Scatter(
            x=sub['N_reads'], y=sub['Expected_genera'],
            mode='lines+markers', name=sample,
            line=dict(color=color, width=2), marker=dict(size=5),
            hovertemplate=(
                f'<b>{sample}</b><br>'
                'reads = %{x:,}<br>'
                'E[genera] = %{y:.2f}<extra></extra>')))
    fig.update_layout(
        title=dict(text='Rarefaction curves · expected unique genera vs '
                        'subsample size (Bracken, Hurlbert closed form)',
                   font_color=BLUE, font_size=14),
        xaxis=dict(title='Subsample size (reads)', type='log'),
        yaxis=dict(title='Expected number of distinct genera'),
        height=380, margin=dict(l=70, r=30, t=60, b=60),
        legend=dict(font_size=9),
        font=dict(family='Arial, sans-serif', size=11))
    return fig


def chart_whittaker(bracken_df, top_n=50):
    """Rank-abundance (Whittaker) plot: log10(abundance) vs rank, one
    line per sample. The slope tells you whether the community is
    dominated by a few taxa (steep) or balanced (shallow)."""
    if bracken_df is None or bracken_df.empty:
        return None
    fig = go.Figure()
    samples = sorted(bracken_df['Sample'].unique())
    for i, sample in enumerate(samples):
        sub = bracken_df[bracken_df['Sample'] == sample].copy()
        sub = sub[sub['Pct'] > 0].sort_values('Pct', ascending=False).head(top_n)
        if sub.empty:
            continue
        sub['Rank'] = np.arange(1, len(sub) + 1)
        color = PALETTE[i % len(PALETTE)]
        fig.add_trace(go.Scatter(
            x=sub['Rank'], y=sub['Pct'],
            mode='lines+markers', name=sample,
            line=dict(color=color, width=2), marker=dict(size=6),
            text=sub['Genus'],
            hovertemplate=(
                f'<b>{sample}</b><br>'
                'rank %{x}: <b>%{text}</b><br>'
                '%{y:.4f}%<extra></extra>')))
    fig.update_layout(
        title=dict(text=f'Whittaker rank-abundance · top-{top_n} genera per sample',
                   font_color=BLUE, font_size=14),
        xaxis=dict(title='Rank (most → least abundant)'),
        yaxis=dict(title='Relative abundance (%)', type='log'),
        height=380, margin=dict(l=70, r=30, t=60, b=60),
        legend=dict(font_size=9),
        font=dict(family='Arial, sans-serif', size=11))
    return fig


def load_phenotypic_validation(results_dir):
    """Load the per-sample phenotypic validation TSVs produced by
    minimap2 against the surveillance target genomes. Columns per file:
    sample · species · genome_size · mapped · mean_depth ·
    breadth_1x · breadth_10x · breadth_30x."""
    rows = []
    for f in sorted(glob.glob(
            f'{results_dir}/12_phenotypic_minimap2/*_phenotypic.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
        except Exception:
            continue
        if df.empty:
            continue
        for _, r in df.iterrows():
            rows.append({
                'Sample'      : str(r.get('sample', '')).strip(),
                'Species'     : str(r.get('species', '')).strip(),
                'Genome_size' : int(r.get('genome_size', 0) or 0),
                'Mapped'      : int(r.get('mapped', 0) or 0),
                'Mean_depth'  : float(r.get('mean_depth', 0) or 0),
                'Breadth_1x'  : float(r.get('breadth_1x', 0) or 0),
                'Breadth_10x' : float(r.get('breadth_10x', 0) or 0),
                'Breadth_30x' : float(r.get('breadth_30x', 0) or 0),
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def html_phenotypic_table(pheno_df):
    """Triage table for whole-genome phenotypic verification.

    The verdict column is a report-side inference: ANY species with
    breadth_10x ≥ 30% AND mean_depth ≥ 5× is treated as "Confirmed
    presence" — both signals must align, because high breadth at 1× is
    just chance contact while low breadth at high depth is normally a
    repeated region. Numbers are conservative; all underlying values are
    in the table."""
    if pheno_df is None or pheno_df.empty:
        return ('<p class="no-data">No phenotypic minimap2 output for '
                'this run (12_phenotypic_minimap2/*_phenotypic.tsv missing). '
                'Disable with <code>--phenotypic_targets false</code> if '
                'this is intended.</p>')
    df = pheno_df.copy()
    df['_priority'] = -df['Mean_depth'] * df['Breadth_10x']
    df = df.sort_values(['Sample', '_priority'])
    rows = ''
    for _, r in df.iterrows():
        depth = float(r['Mean_depth'])
        b10 = float(r['Breadth_10x'])
        b30 = float(r['Breadth_30x'])
        if depth >= 5.0 and b10 >= 30.0:
            verdict_color, verdict_label = GREEN, 'Confirmed'
        elif depth >= 1.0 and b10 >= 5.0:
            verdict_color, verdict_label = ORANGE, 'Partial'
        elif r['Mapped'] > 0:
            verdict_color, verdict_label = LBLUE, 'Trace mapping'
        else:
            verdict_color, verdict_label = GREY, 'Absent'
        depth_cls = ('good' if depth >= 10 else
                     ('warn' if depth >= 1 else 'note'))
        b10_cls = ('good' if b10 >= 50 else
                   ('warn' if b10 >= 10 else 'note'))
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["Sample"])}">'
            f'{_esc(r["Sample"])}</td>'
            f'<td><em>{_esc(r["Species"])}</em></td>'
            f'<td>{r["Genome_size"]:,}</td>'
            f'<td>{r["Mapped"]:,}</td>'
            f'<td class="{depth_cls}">{depth:.2f}×</td>'
            f'<td>{r["Breadth_1x"]:.2f}%</td>'
            f'<td class="{b10_cls}">{b10:.2f}%</td>'
            f'<td>{b30:.2f}%</td>'
            f'<td><span class="badge" style="background:{verdict_color}">'
            f'{verdict_label}</span></td></tr>')
    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Target species</th>'
        '<th>Genome size (bp)</th><th>Reads mapped</th>'
        '<th>Mean depth</th><th>Breadth ≥1×</th>'
        '<th>Breadth ≥10×</th><th>Breadth ≥30×</th>'
        '<th>Verdict</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note"><strong>From the data</strong>: depth and breadth '
        'columns come straight from <code>samtools coverage</code> on the '
        'phenotypic-targets minimap2 alignment. '
        '<strong>Inferred by this report</strong>: the Verdict badge — '
        '<em>Confirmed</em> = mean depth ≥ 5× AND breadth ≥10× ≥ 30% '
        '(both signals aligned, real organism); '
        '<em>Partial</em> = mean depth ≥ 1× AND breadth ≥10× ≥ 5% '
        '(probable, low coverage); '
        '<em>Trace mapping</em> = some reads mapped but neither threshold '
        'met (chance hits, related organism, or homologous region); '
        '<em>Absent</em> = zero mapped reads. '
        'Rows sorted within each sample by depth × breadth descending.</p>')


def load_amrfinder_catalog(amrfinder_db_path):
    """Load the AMRFinderPlus reference-gene catalogue (fam.tsv) into a
    lookup keyed by NF identifier (with and without the `.N` version
    suffix, so callers can query either form).

    Each entry: {gene_symbol, family, product, type, subtype, class, subclass}.

    fam.tsv columns (AMRFinderPlus 4.x):
      node_id · parent_node_id · gene_symbol · hmm_id · …
      type · subtype · class · subclass · family_name

    The NF identifier that geNomad/AMRFinderPlus print externally lives
    in the `hmm_id` column (`NF033135.1`); we strip the version when
    indexing so `NF033135` and `NF033135.1` both resolve.
    """
    if not amrfinder_db_path:
        return {}
    candidates = [
        os.path.join(amrfinder_db_path, 'fam.tsv'),
        # Older AMRFinderPlus layouts kept it one level up
        os.path.join(amrfinder_db_path, '..', 'fam.tsv'),
    ]
    fam_path = next((p for p in candidates if os.path.isfile(p)), None)
    if fam_path is None:
        return {}
    catalog = {}
    try:
        with open(fam_path) as fh:
            header = fh.readline().lstrip('#').rstrip('\n').split('\t')
            cols = {name: i for i, name in enumerate(header)}
        df = pd.read_csv(fam_path, sep='\t', comment='#', header=None,
                         names=header, keep_default_na=False, dtype=str)
        for _, r in df.iterrows():
            hmm_id = str(r.get('hmm_id', '')).strip()
            if not hmm_id or hmm_id == '-':
                continue
            base = hmm_id.split('.', 1)[0]
            entry = {
                'nf_id'       : hmm_id,
                'gene_symbol' : str(r.get('gene_symbol', '') or '').strip(),
                'type'        : str(r.get('type', '') or '').strip(),
                'subtype'     : str(r.get('subtype', '') or '').strip(),
                'class'       : str(r.get('class', '') or '').strip(),
                'subclass'    : str(r.get('subclass', '') or '').strip(),
                'family'      : str(r.get('family_name', '') or '').strip(),
                'node_id'     : str(r.get('node_id', '') or '').strip(),
            }
            catalog.setdefault(base, entry)
            catalog.setdefault(hmm_id, entry)
    except Exception as e:
        print(f'[WARN] Could not parse {fam_path}: {e}')
    return catalog


def load_run_params(results_dir):
    """Load the params.json emitted by `workflow.onComplete`. Returns
    {} if absent so the report still renders for older runs."""
    import json
    candidates = [
        f'{results_dir}/00_reports/params.json',
        f'{results_dir}/params.json',
    ]
    for p in candidates:
        if os.path.isfile(p):
            try:
                with open(p) as fh:
                    return json.load(fh)
            except Exception as e:
                print(f'[WARN] Could not parse {p}: {e}')
    return {}


# ═══════════════════════════════════════════════════════════════
# BRANDING
# ═══════════════════════════════════════════════════════════════

def load_branding(branding_dir):
    info = {'name': 'EpiTaxMAG', 'department': '', 'address': '',
            'contact': '', 'phone': '', 'web': ''}
    logo_b64 = ''

    info_path = os.path.join(branding_dir, 'org.info')
    if os.path.isfile(info_path):
        with open(info_path) as f:
            for line in f:
                if '=' in line:
                    k, v = line.strip().split('=', 1)
                    info[k.strip()] = v.strip()

    for ext in ['png', 'svg', 'jpg', 'ico']:
        logo_path = os.path.join(branding_dir, f'logo.{ext}')
        if os.path.isfile(logo_path):
            mime = {'svg': 'image/svg+xml', 'png': 'image/png',
                    'ico': 'image/x-icon', 'jpg': 'image/jpeg'}[ext]
            with open(logo_path, 'rb') as f:
                logo_b64 = f'data:{mime};base64,{base64.b64encode(f.read()).decode()}'
            break

    return info, logo_b64


# ═══════════════════════════════════════════════════════════════
# DATA LOADERS
# ═══════════════════════════════════════════════════════════════

def parse_float(val):
    if isinstance(val, (int, float)):
        return float(val)
    return float(str(val).replace(',', '.'))


def load_qc_csv(path):
    if not os.path.isfile(path):
        return pd.DataFrame()
    with open(path) as f:
        h = f.readline()
    sep = '\t' if '\t' in h else (';' if ';' in h else ',')
    df = pd.read_csv(path, sep=sep)
    for c in df.columns:
        if c != 'Sample':
            df[c] = df[c].apply(lambda v: parse_float(v) if pd.notna(v) else 0)
    return df


def load_run_stats(input_dir):
    """Count reads from input fastq.gz files including unclassified."""
    import gzip
    stats = {'samples': [], 'total_reads': 0, 'classified_reads': 0,
             'unclassified_reads': 0, 'pct_unclassified': 0}
    if not input_dir or not os.path.isdir(input_dir):
        return stats
    for f in sorted(glob.glob(f'{input_dir}/*.fastq.gz')):
        name = os.path.basename(f).replace('.fastq.gz', '')
        try:
            n = 0
            with gzip.open(f, 'rt') as fh:
                for line in fh:
                    if line.startswith('@'):
                        n += 1
            stats['total_reads'] += n
            if 'unclassified' in name.lower():
                stats['unclassified_reads'] = n
            else:
                stats['classified_reads'] += n
                stats['samples'].append({'name': name, 'reads': n})
        except Exception:
            pass
    if stats['total_reads'] > 0:
        stats['pct_unclassified'] = round(
            stats['unclassified_reads'] / stats['total_reads'] * 100, 1)
    return stats


def load_checkm2(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/23_binqc_checkm2/*/quality_report.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                comp = float(r['Completeness'])
                cont = float(r['Contamination'])
                qc = 'HQ' if comp >= 90 and cont <= 5 else (
                    'MQ' if comp >= 50 and cont <= 10 else 'LQ')
                rows.append({
                    'Sample': sample, 'Bin': r['Name'],
                    'Completeness': comp, 'Contamination': cont,
                    'Quality': qc,
                    'Genome_Size': int(r.get('Genome_Size', 0)),
                    'GC_Content': float(r.get('GC_Content', 0)),
                    'Contig_N50': int(r.get('Contig_N50', 0)),
                    'Total_Contigs': int(r.get('Total_Contigs', 0)),
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


TAXONOMY_TOOLS = [
    # (subdir, label, version, db_label)
    ('24_taxonomy_skani',    'Skani',    '0.3.1', 'GTDB r226'),
    ('24_taxonomy_sourmash', 'Sourmash', '4.9.5', 'GTDB r226'),
    ('24_taxonomy_gtdbtk',   'GTDB-Tk',  '2.7.0', 'GTDB r232'),
]


def detect_taxonomy_tool(results_dir):
    """Detect which taxonomy classifier produced the 24_taxonomy_* output.
    Returns dict with subdir, label, version, db_label, full label/short label."""
    for subdir, label, version, db_label in TAXONOMY_TOOLS:
        if glob.glob(f'{results_dir}/{subdir}/*/*_taxonomy.tsv'):
            return {
                'subdir': subdir,
                'label': label,
                'version': version,
                'db_label': db_label,
                'full': f'{label} {version} ({db_label})',
                'short': f'{label} ({db_label})',
            }
    # Default fallback for empty/missing taxonomy
    return {
        'subdir': '24_taxonomy_skani',
        'label': 'Skani',
        'version': '0.3.1',
        'db_label': 'GTDB r226',
        'full': 'Skani 0.3.1 (GTDB r226)',
        'short': 'Skani (GTDB r226)',
    }


def _split_skani_classification(cls):
    """Skani classification column: '<accession> <genus species ... description>'.
    Try to extract genus/species from the free text."""
    cls = cls.strip()
    if not cls:
        return '', '', '', ''
    parts = cls.split(None, 1)
    accession = parts[0] if len(parts) > 0 else ''
    rest = parts[1] if len(parts) > 1 else ''
    # Stop tokens that mark the end of the binomial in NCBI titles
    stop_tokens = (' strain ', ' isolate ', ' sp. ', ' subsp. ', ' DSM ',
                   ' ATCC ', ' CIP ', ' NCTC ', ' contig', ' scaffold',
                   ' chromosome', ' complete', ' whole genome', ' MAG:', ',')
    cut = len(rest)
    for tok in stop_tokens:
        i = rest.find(tok)
        if i >= 0 and i < cut:
            cut = i
    binomial = rest[:cut].strip()
    tokens = binomial.split()
    genus = tokens[0] if tokens else ''
    species = ' '.join(tokens[1:]) if len(tokens) > 1 else ''
    return accession, genus, species, binomial


def _split_gtdb_classification(cls):
    """GTDB-Tk style 'd__...;p__...;c__...;o__...;f__...;g__...;s__...'."""
    levels = {}
    for part in cls.split(';'):
        p = part.strip()
        for px, lv in [('d__', 'domain'), ('p__', 'phylum'), ('c__', 'class'),
                       ('o__', 'order'), ('f__', 'family'), ('g__', 'genus'),
                       ('s__', 'species')]:
            if p.startswith(px):
                levels[lv] = p[len(px):] or ''
    return levels


def load_taxonomy(results_dir, tax_info=None):
    if tax_info is None:
        tax_info = detect_taxonomy_tool(results_dir)
    subdir = tax_info['subdir']
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/{subdir}/*/*_taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue
            for _, r in df.iterrows():
                bin_name = str(r.iloc[0]).strip()
                cls = str(r.iloc[1]) if len(r) > 1 else ''
                if subdir == '24_taxonomy_gtdbtk':
                    levels = _split_gtdb_classification(cls)
                    org = (f"{levels.get('genus', '')} "
                           f"{levels.get('species', '')}").strip() or 'Unclassified'
                    rows.append({
                        'Sample': sample, 'Bin': bin_name,
                        'Taxonomy': cls, 'Organism': org,
                        'Domain': levels.get('domain', ''),
                        'Phylum': levels.get('phylum', ''),
                        'Class': levels.get('class', ''),
                        'Order': levels.get('order', ''),
                        'Family': levels.get('family', ''),
                        'Genus': levels.get('genus', ''),
                        'Species': levels.get('species', ''),
                        'Reference': '', 'ANI': np.nan, 'AF': np.nan,
                    })
                else:
                    # Skani / Sourmash: free-text classification with leading accession
                    accession, genus, species, binomial = _split_skani_classification(cls)
                    org = binomial if binomial else 'Unclassified'
                    try:
                        ani = float(r.get('ani', np.nan))
                    except (ValueError, TypeError):
                        ani = np.nan
                    try:
                        af = float(r.get('af', np.nan))
                    except (ValueError, TypeError):
                        af = np.nan
                    rows.append({
                        'Sample': sample, 'Bin': bin_name,
                        'Taxonomy': cls, 'Organism': org,
                        'Domain': '', 'Phylum': '',
                        # Use genus as 'Class' fallback so the colored stacked
                        # chart still groups meaningfully (Skani/Sourmash do
                        # not return a class-rank string).
                        'Class': genus, 'Order': '', 'Family': '',
                        'Genus': genus, 'Species': species,
                        'Reference': accession, 'ANI': ani, 'AF': af,
                    })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_amr_mags(results_dir):
    rows = []
    for per_bin_dir in sorted(glob.glob(f'{results_dir}/25_amr_mags/*/per_bin')):
        sample = os.path.basename(os.path.dirname(per_bin_dir))
        for f in sorted(glob.glob(f'{per_bin_dir}/*_amr.tsv')):
            try:
                df = pd.read_csv(f, sep='\t')
                for _, r in df.iterrows():
                    rows.append({
                        'Sample': sample,
                        'Bin': str(r.get('Name', '')).strip(),
                        'Contig': str(r.get('Contig id', '')).strip(),
                        'Gene': str(r.get('Element symbol', r.get('Gene symbol', ''))).strip(),
                        'Gene_desc': str(r.get('Element name', r.get('Sequence name', ''))).strip(),
                        'Scope': str(r.get('Scope', '')).strip(),
                        'Type': str(r.get('Type', r.get('Element type', ''))).strip(),
                        'Subtype': str(r.get('Subtype', r.get('Element subtype', ''))).strip(),
                        'Class': str(r.get('Class', '')).strip(),
                        'Subclass': str(r.get('Subclass', '')).strip(),
                        'Method': str(r.get('Method', '')).strip(),
                        'Identity': float(r.get('% Identity to reference',
                                                 r.get('% Identity to reference sequence', 0))),
                        'Coverage': float(r.get('% Coverage of reference',
                                                r.get('% Coverage of reference sequence', 0))),
                        'Accession': str(r.get('Closest reference accession', '')).strip(),
                        'Closest_ref': str(r.get('Closest reference name', '')).strip(),
                        'HMM_desc': str(r.get('HMM description', '')).strip(),
                    })
            except Exception:
                pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_kma(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/11_amr_kma/*.res')):
        sample = os.path.basename(f).replace('.res', '')
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue
            for _, r in df.iterrows():
                template = str(r.iloc[0]).strip()
                if not template or template.startswith('#'):
                    continue
                gene_root = template.split('_')[0] if '_' in template else template
                rows.append({
                    'Sample': sample, 'Template': template, 'Gene_root': gene_root,
                    'KMA_Identity': float(r.get('Template_Identity', 0)),
                    'KMA_Coverage': float(r.get('Template_Coverage', 0)),
                    'KMA_Depth': float(r.get('Depth', 0)),
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_plasmids(results_dir):
    """geNomad plasmid_summary loader. Note: pandas.read_csv treats the
    literal string `NA` from geNomad as NaN by default — we keep it
    literal here so downstream "is this empty?" checks work without
    isnull contortions."""
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/26_genomad/*/*_plasmid_summary.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t', keep_default_na=False, na_values=[''])
            for _, r in df.iterrows():
                contig = str(r.iloc[0]).strip().split('|')[0]
                score = float(r.get('plasmid_score', r.get('score', 0)) or 0)
                n_genes = int(r.get('n_genes', 0) or 0)
                amr_raw = r.get('amr_genes', '')
                amr_genes = '' if str(amr_raw).strip() in ('', 'NA', 'nan') \
                    else str(amr_raw).strip()
                rows.append({
                    'Sample': sample, 'Contig': contig,
                    'Plasmid_score': score,
                    'Length': int(r.get('length', 0) or 0),
                    'N_genes': n_genes,
                    'AMR_genes': amr_genes,
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_amr_plasmids(results_dir):
    """AMRFinderPlus on geNomad plasmid contigs (Phase B output)."""
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/26_genomad/*/*_plasmid_amr.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue
            for _, r in df.iterrows():
                rows.append({
                    'Sample': str(r.get('Sample', '')).strip(),
                    'Contig': str(r.get('Contig', '')).strip(),
                    'Gene': str(r.get('Gene', '')).strip(),
                    'Name': str(r.get('Name', '')).strip(),
                    'Class': str(r.get('Class', '')).strip(),
                    'Subclass': str(r.get('Subclass', '')).strip(),
                    'Method': str(r.get('Method', '')).strip(),
                    'Identity': float(r.get('Identity', 0) or 0),
                    'Coverage': float(r.get('Coverage', 0) or 0),
                    'Accession': str(r.get('Accession', '')).strip(),
                    'Closest_ref': str(r.get('Closest_ref', '')).strip(),
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_integration(results_dir):
    dfs = []
    for f in sorted(glob.glob(f'{results_dir}/27_amr_pathogen_integration/*_amr_pathogen.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty:
                dfs.append(df)
        except Exception:
            pass
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def load_vfdb(results_dir):
    """Load VFDB ABRicate results. Header starts with #FILE."""
    rows = []
    VFDB_COLUMNS = ['FILE', 'SEQUENCE', 'START', 'END', 'STRAND', 'GENE',
                     'COVERAGE', 'COVERAGE_MAP', 'GAPS', 'PCT_COVERAGE',
                     'PCT_IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT',
                     'RESISTANCE']
    for f in sorted(glob.glob(f'{results_dir}/28_abricate/*/*.vfdb.tsv')):
        sample_dir = os.path.basename(os.path.dirname(f))
        bin_name = os.path.basename(f).replace('.vfdb.tsv', '')
        # Remove sample prefix from bin name if present
        if bin_name.startswith(sample_dir + '_'):
            bin_short = bin_name[len(sample_dir) + 1:]
        else:
            bin_short = bin_name
        try:
            df = pd.read_csv(f, sep='\t', header=0)
            df.columns = VFDB_COLUMNS[:len(df.columns)]
            # Skip header row if first column value starts with '#'
            df = df[~df['FILE'].astype(str).str.startswith('#')]
            if df.empty:
                continue
            for _, r in df.iterrows():
                identity = 0
                coverage = 0
                try:
                    identity = float(r.get('PCT_IDENTITY', 0))
                except (ValueError, TypeError):
                    pass
                try:
                    coverage = float(r.get('PCT_COVERAGE', 0))
                except (ValueError, TypeError):
                    pass
                rows.append({
                    'Sample': sample_dir,
                    'Bin': bin_short,
                    'Gene': str(r.get('GENE', '')),
                    'Product': str(r.get('PRODUCT', '')),
                    'Accession': str(r.get('ACCESSION', '')),
                    'Identity': identity,
                    'Coverage': coverage,
                    'Sequence': str(r.get('SEQUENCE', '')),
                    'Database': str(r.get('DATABASE', '')),
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_quast(results_dir):
    """Load QUAST report.tsv (transposed: samples are columns)."""
    path = f'{results_dir}/17_qc_assembly/quast_results/report.tsv'
    if not os.path.isfile(path):
        return pd.DataFrame()
    try:
        raw = pd.read_csv(path, sep='\t', index_col=0)
        # Transpose so rows are samples
        df = raw.T.reset_index()
        df.rename(columns={'index': 'Sample'}, inplace=True)
        return df
    except Exception:
        return pd.DataFrame()


def count_bins(results_dir):
    """Count bins per binner per sample."""
    binners = {
        'MetaBAT2': f'{results_dir}/19_binning_metabat2/*/metabat_bins/*.fa',
        'MaxBin2': f'{results_dir}/20_binning_maxbin2/*/maxbin_bins/*.fa',
        'SemiBin2': f'{results_dir}/21_binning_semibin2/*/semibin_bins/*.fa',
        'DAS Tool': f'{results_dir}/22_binning_dastool/*/*_DASTool_bins/*.fa',
    }
    rows = []
    for binner, pattern in binners.items():
        files = glob.glob(pattern)
        per_sample = {}
        for f in files:
            # Extract sample name from path
            parts = f.split('/')
            for i, p in enumerate(parts):
                if p in ('19_binning_metabat2', '20_binning_maxbin2',
                         '21_binning_semibin2', '22_binning_dastool'):
                    sample = parts[i + 1] if i + 1 < len(parts) else ''
                    break
            else:
                continue
            per_sample[sample] = per_sample.get(sample, 0) + 1
        for sample, count in sorted(per_sample.items()):
            rows.append({'Sample': sample, 'Binner': binner, 'Bins': count})
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_bakta_summaries(results_dir):
    """Parse Bakta .txt summary files for tRNA/rRNA/CDS counts per MAG."""
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/31_bakta/*/*/*.txt')):
        parts = f.split('/')
        # Structure: .../31_bakta/SAMPLE/BIN/BIN.txt
        bin_name = os.path.basename(os.path.dirname(f))
        sample = parts[-3] if len(parts) >= 3 else ''
        try:
            data = {'Sample': sample, 'Bin': bin_name}
            with open(f) as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith('Length:'):
                        data['Length'] = int(line.split(':')[1].strip())
                    elif line.startswith('Count:'):
                        data['Seq_count'] = int(line.split(':')[1].strip())
                    elif line.startswith('GC:'):
                        data['GC'] = float(line.split(':')[1].strip())
                    elif line.startswith('N50:'):
                        data['N50'] = int(line.split(':')[1].strip())
                    elif line.startswith('tRNAs:'):
                        data['tRNAs'] = int(line.split(':')[1].strip())
                    elif line.startswith('tmRNAs:'):
                        data['tmRNAs'] = int(line.split(':')[1].strip())
                    elif line.startswith('rRNAs:'):
                        data['rRNAs'] = int(line.split(':')[1].strip())
                    elif line.startswith('ncRNAs:'):
                        data['ncRNAs'] = int(line.split(':')[1].strip())
                    elif line.startswith('CRISPR arrays:'):
                        data['CRISPRs'] = int(line.split(':')[1].strip())
                    elif line.startswith('CDSs:'):
                        data['CDSs'] = int(line.split(':')[1].strip())
                    elif line.startswith('pseudogenes:'):
                        data['Pseudogenes'] = int(line.split(':')[1].strip())
                    elif line.startswith('hypotheticals:'):
                        data['Hypotheticals'] = int(line.split(':')[1].strip())
                    elif line.startswith('coding density:'):
                        try:
                            data['Coding_density'] = float(line.split(':')[1].strip())
                        except ValueError:
                            pass
                    elif line.startswith('oriCs:'):
                        data['oriCs'] = int(line.split(':')[1].strip())
                    elif line.startswith('oriTs:'):
                        data['oriTs'] = int(line.split(':')[1].strip())
            rows.append(data)
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_integrons(results_dir):
    """Load IntegronFinder results. Filter to meaningful rows only."""
    rows = []
    for f in sorted(glob.glob(
            f'{results_dir}/30_integronfinder/*/*/Results_*/*.integrons')):
        # Derive sample and bin from path
        # .../30_integronfinder/SAMPLE/BIN/Results_*/FILE.integrons
        parts = f.split('/')
        idx = None
        for i, p in enumerate(parts):
            if p == '30_integronfinder':
                idx = i
                break
        if idx is None:
            continue
        sample = parts[idx + 1] if idx + 1 < len(parts) else ''
        bin_name = parts[idx + 2] if idx + 2 < len(parts) else ''

        try:
            with open(f) as fh:
                lines = fh.readlines()
            header = None
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    # Check if this is a "No Integron found" message
                    if 'No Integron' in line:
                        break
                    continue
                if header is None:
                    # First non-comment non-empty line should be header
                    # But the header line starts with ID_integron
                    if line.startswith('ID_integron'):
                        header = line.split('\t')
                        continue
                    else:
                        # Try to detect header from the file
                        # Parse anyway
                        header = ['ID_integron', 'ID_replicon', 'element',
                                  'pos_beg', 'pos_end', 'strand', 'evalue',
                                  'type_elt', 'annotation', 'model', 'type',
                                  'default', 'distance_2attC',
                                  'considered_topology']
                fields = line.split('\t')
                if len(fields) < 11:
                    continue
                row = dict(zip(header, fields))
                row['Sample'] = sample
                row['Bin'] = bin_name
                rows.append(row)
        except Exception:
            pass

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    return df


def load_mobsuite(results_dir):
    """Load MOB-suite mobtyper results."""
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/29_mobsuite/*/mobtyper_results.txt')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample,
                    'sample_id': str(r.get('sample_id', '')),
                    'num_contigs': int(r.get('num_contigs', 0)),
                    'size': int(r.get('size', 0)),
                    'gc': float(r.get('gc', 0)),
                    'rep_type': str(r.get('rep_type(s)', '')),
                    'relaxase_type': str(r.get('relaxase_type(s)', '')),
                    'mpf_type': str(r.get('mpf_type', '')),
                    'predicted_mobility': str(r.get('predicted_mobility', '')),
                    'mash_nearest_neighbor': str(r.get('mash_nearest_neighbor', '')),
                    'mash_neighbor_identification': str(
                        r.get('mash_neighbor_identification', '')),
                    'primary_cluster_id': str(r.get('primary_cluster_id', '')),
                })
        except Exception:
            pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ═══════════════════════════════════════════════════════════════
# ANALYSIS
# ═══════════════════════════════════════════════════════════════

def compute_survival(qc_raw, qc_filt):
    if qc_raw.empty or qc_filt.empty:
        return pd.DataFrame()
    raw = qc_raw[['Sample', 'total_reads', 'mean_length']].copy()
    raw.columns = ['Sample', 'raw_reads', 'raw_mean_len']
    raw['raw_gb'] = raw['raw_reads'] * raw['raw_mean_len'] / 1e9
    filt = qc_filt[['Sample', 'total_reads', 'mean_length']].copy()
    filt.columns = ['Sample', 'clean_reads', 'clean_mean_len']
    filt['clean_gb'] = filt['clean_reads'] * filt['clean_mean_len'] / 1e9
    df = pd.merge(raw, filt, on='Sample', how='outer').fillna(0)
    df['pct_reads'] = np.where(
        df['raw_reads'] > 0, df['clean_reads'] / df['raw_reads'] * 100, 0
    ).round(1)
    df['pct_gb'] = np.where(
        df['raw_gb'] > 0, df['clean_gb'] / df['raw_gb'] * 100, 0
    ).round(1)
    return df.sort_values('Sample').reset_index(drop=True)


def compute_concordance(kma_df, amr_df):
    if kma_df.empty and amr_df.empty:
        return pd.DataFrame()
    rows = []
    samples = set()
    if not kma_df.empty:
        samples.update(kma_df['Sample'].unique())
    if not amr_df.empty:
        samples.update(amr_df['Sample'].unique())

    for sample in sorted(samples):
        kma_genes = set()
        kma_info = {}
        if not kma_df.empty:
            sk = kma_df[kma_df['Sample'] == sample]
            for _, r in sk.iterrows():
                root = r['Gene_root']
                kma_genes.add(root)
                if root not in kma_info or r['KMA_Depth'] > kma_info[root]['KMA_Depth']:
                    kma_info[root] = {
                        'KMA_Depth': r['KMA_Depth'],
                        'KMA_Identity': r['KMA_Identity'],
                        'KMA_Coverage': r['KMA_Coverage'],
                    }

        mag_genes = set()
        mag_info = {}
        if not amr_df.empty:
            sa = amr_df[amr_df['Sample'] == sample]
            for _, r in sa.iterrows():
                gene = r['Gene']
                mag_genes.add(gene)
                if gene not in mag_info:
                    mag_info[gene] = {
                        'MAG_Identity': r['Identity'],
                        'MAG_Coverage': r['Coverage'],
                        'Bin': r['Bin'],
                        'Class': r['Class'],
                    }

        all_genes = kma_genes | mag_genes
        for gene in sorted(all_genes):
            in_kma = gene in kma_genes
            in_mag = gene in mag_genes
            source = 'BOTH' if in_kma and in_mag else (
                'READS' if in_kma else 'CONTIGS')
            row = {'Sample': sample, 'Gene': gene, 'Source': source}
            if in_kma:
                row.update(kma_info.get(gene, {}))
            if in_mag:
                row.update(mag_info.get(gene, {}))
            rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def html_risk_model_panel():
    """Render the risk-model decision logic verbatim from
    integrate_amr_pathogen.RISK_RULES, plus the high-risk class list.
    Imported lazily so a missing module never blocks the report."""
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from integrate_amr_pathogen import RISK_RULES, RISK_HIGH_CLASSES
    except Exception as e:
        return f'<p class="no-data">Risk model not available: {e}</p>'

    rule_rows = ''
    badge_html = {
        'CRITICAL': RISK_BADGE['CRITICAL'],
        'HIGH'    : RISK_BADGE['HIGH'],
        'MEDIUM'  : RISK_BADGE['MEDIUM'],
        'LOW'     : RISK_BADGE['LOW'],
    }
    for rule_id, reason, _pred, level in RISK_RULES:
        rule_rows += (
            f'<tr><td><strong>{_esc(rule_id)}</strong></td>'
            f'<td>{badge_html.get(level, "")}</td>'
            f'<td>{_esc(reason)}</td></tr>')

    classes_html = ', '.join(
        f'<code>{_esc(c)}</code>' for c in sorted(RISK_HIGH_CLASSES))

    return (
        '<div class="risk-model-panel">'
        '<p>Risk evaluation is a deterministic, ordered set of rules; '
        'the first rule whose predicate matches a given AMR call '
        'sets the level. Each row in the AMR detail table carries the '
        'matched rule ID so the decision is fully traceable.</p>'
        '<h4>High-risk antibiotic classes</h4>'
        f'<p class="risk-classes">{classes_html}</p>'
        '<p class="note">Membership of this set is what flips a call '
        'between LOW/MEDIUM and HIGH/CRITICAL when the location is '
        'a confirmed plasmid or the MAG is HQ.</p>'
        '<h4>Rules (evaluated top-to-bottom)</h4>'
        '<table class="prov-table"><thead><tr>'
        '<th style="width:60px">Rule</th>'
        '<th style="width:100px">Level</th>'
        '<th>Trigger</th></tr></thead>'
        f'<tbody>{rule_rows}</tbody></table>'
        '</div>')


def confidence_score(identity, coverage, method, kma_concordance,
                     mag_quality, plasmid_score):
    score = 0
    score += min(identity, 100) * 0.3
    score += min(coverage, 100) * 0.2
    method_scores = {'EXACT': 20, 'EXACTX': 20, 'ALLELE': 18, 'BLAST': 15,
                     'PARTIAL': 10, 'HMM': 12, 'INTERNAL_STOP': 8,
                     'PARTIAL_CONTIG_END': 10, 'POINT': 18}
    score += method_scores.get(method, 10)
    if kma_concordance:
        score += 15
    qc_scores = {'HQ': 15, 'MQ': 10, 'LQ': 5}
    score += qc_scores.get(mag_quality, 0)
    return min(round(score), 100)


def summarize_integrons(integron_df):
    """Summarize integrons per sample/bin. Filter to integrons that contain intI.

    Only keep rows belonging to an ID_integron group that has at least one
    intI annotation, then count types (complete, In0, CALIN) from those
    filtered integrons.
    """
    if integron_df.empty:
        return pd.DataFrame()

    # Step 1: Identify integron IDs that contain at least one intI
    has_annotation = 'annotation' in integron_df.columns
    has_id = 'ID_integron' in integron_df.columns
    has_type = 'type' in integron_df.columns

    if has_annotation and has_id:
        intI_rows = integron_df[integron_df['annotation'] == 'intI']
        # Get the set of (Sample, Bin, ID_integron) that have intI
        intI_keys = set(zip(intI_rows['Sample'], intI_rows['Bin'],
                            intI_rows['ID_integron']))
        # Filter to rows from those integrons only
        df = integron_df[integron_df.apply(
            lambda r: (r['Sample'], r['Bin'], r['ID_integron']) in intI_keys,
            axis=1)].copy()
    else:
        df = integron_df.copy()

    if df.empty:
        return pd.DataFrame()

    rows = []
    groups = df.groupby(['Sample', 'Bin'])
    for (sample, bin_name), grp in groups:
        n_complete = int((grp['type'] == 'complete').sum()) if has_type else 0
        n_in0 = int((grp['type'] == 'In0').sum()) if has_type else 0
        n_calin = int((grp['type'] == 'CALIN').sum()) if has_type else 0
        has_intI = 'Yes' if has_annotation and (
            grp['annotation'] == 'intI').any() else 'No'
        # Only include rows where there is something interesting
        if n_complete > 0 or n_in0 > 0 or n_calin > 0 or has_intI == 'Yes':
            rows.append({
                'Sample': sample,
                'Bin': bin_name,
                'N_complete': n_complete,
                'N_In0': n_in0,
                'N_CALIN': n_calin,
                'has_intI': has_intI,
                'Total_elements': len(grp),
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ═══════════════════════════════════════════════════════════════
# HELPER
# ═══════════════════════════════════════════════════════════════

def _safe_str(val):
    if val is None:
        return ''
    if isinstance(val, float) and np.isnan(val):
        return ''
    s = str(val).strip()
    return '' if s.lower() in ('nan', 'na', 'none', '') else s


def _esc(val):
    """HTML-escape a string."""
    s = _safe_str(val)
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


# ═══════════════════════════════════════════════════════════════
# PLOTLY CHARTS
# ═══════════════════════════════════════════════════════════════

def fig_to_html(fig):
    if fig is None:
        return '<p class="no-data">Sin datos disponibles.</p>'
    return fig.to_html(full_html=False, include_plotlyjs=False)


def chart_retention(surv):
    if surv.empty:
        return None
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=['Reads', 'Gigabases'],
                        horizontal_spacing=0.12)
    s = surv['Sample'].values
    fig.add_trace(go.Bar(
        name='Raw', x=s, y=surv['raw_reads'], marker_color=BLUE, opacity=0.7,
        text=[f'{v:,.0f}' for v in surv['raw_reads']],
        textposition='outside', textfont_size=8), row=1, col=1)
    fig.add_trace(go.Bar(
        name='Clean', x=s, y=surv['clean_reads'], marker_color=GREEN, opacity=0.85,
        text=[f'{v:,.0f}' for v in surv['clean_reads']],
        textposition='outside', textfont_size=8), row=1, col=1)
    fig.add_trace(go.Bar(
        name='Raw Gb', x=s, y=surv['raw_gb'], marker_color=BLUE, opacity=0.7,
        showlegend=False,
        text=[f'{v:.2f}' for v in surv['raw_gb']],
        textposition='outside', textfont_size=8), row=1, col=2)
    fig.add_trace(go.Bar(
        name='Clean Gb', x=s, y=surv['clean_gb'], marker_color=GREEN, opacity=0.85,
        showlegend=False,
        text=[f'{v:.2f}' for v in surv['clean_gb']],
        textposition='outside', textfont_size=8), row=1, col=2)
    fig.update_layout(
        barmode='group', height=400,
        margin=dict(l=60, r=30, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.04,
                    xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif'))
    fig.update_xaxes(tickangle=-40, tickfont_size=8)
    return fig


def chart_summary_bar(data):
    """Summary table + grouped bar chart: per sample counts of MAGs, AMR, VF,
    Plasmids, Integrons (intI count). Uses table + chart to handle scale
    differences between categories."""
    surv = data['survival']
    checkm2 = data['checkm2']
    integ = data['integration']
    vfdb = data['vfdb']
    plasmids = data['plasmids']
    integrons_raw = data.get('integrons_raw', pd.DataFrame())

    if surv.empty:
        return None

    samples = sorted(surv['Sample'].unique())

    mag_counts = []
    amr_counts = []
    vf_counts = []
    plasmid_counts = []
    intI_counts = []

    for s in samples:
        mag_counts.append(
            len(checkm2[checkm2['Sample'] == s]) if not checkm2.empty else 0)
        amr_counts.append(
            len(integ[integ['Sample'] == s]) if not integ.empty else 0)
        vf_counts.append(
            len(vfdb[vfdb['Sample'] == s]) if not vfdb.empty else 0)
        plasmid_counts.append(
            len(plasmids[plasmids['Sample'] == s]) if not plasmids.empty else 0)
        # Use intI (integrase) count only, not raw element count
        if (not integrons_raw.empty and 'annotation' in integrons_raw.columns
                and 'Sample' in integrons_raw.columns):
            s_intI = integrons_raw[
                (integrons_raw['Sample'] == s) &
                (integrons_raw['annotation'] == 'intI')]
            intI_counts.append(len(s_intI))
        else:
            intI_counts.append(0)

    fig = go.Figure()
    fig.add_trace(go.Bar(name='MAGs', x=samples, y=mag_counts,
                         marker_color=BLUE))
    fig.add_trace(go.Bar(name='AMR genes', x=samples, y=amr_counts,
                         marker_color=RED))
    fig.add_trace(go.Bar(name='Virulence factors', x=samples, y=vf_counts,
                         marker_color=ORANGE))
    fig.add_trace(go.Bar(name='Plasmids', x=samples, y=plasmid_counts,
                         marker_color='#9B59B6'))
    fig.add_trace(go.Bar(name='Integrases (intI)', x=samples, y=intI_counts,
                         marker_color='#1ABC9C'))
    fig.update_layout(
        barmode='group', height=400,
        title=dict(text='Per-Sample Summary', font_color=BLUE, font_size=14),
        xaxis_tickangle=-40,
        margin=dict(l=60, r=30, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.04,
                    xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=10))
    return fig


def chart_taxonomy_stacked(tax_df, checkm2_df, tax_info=None):
    """Stacked bar: organisms per sample.
    For GTDB-Tk grouping is by class rank; for Skani/Sourmash (which only
    return a free-text best hit) it falls back to genus.
    """
    if tax_df.empty:
        return None
    # Merge with checkm2 to get only quality-passed MAGs
    df = tax_df.copy()
    if not checkm2_df.empty:
        valid_bins = set(zip(checkm2_df['Sample'], checkm2_df['Bin']))
        df = df[df.apply(lambda r: (r['Sample'], r['Bin']) in valid_bins, axis=1)]
    if df.empty:
        return None

    # For Skani/Sourmash 'Class' was populated with genus, so this still works.
    df['Color_group'] = df['Class'].apply(
        lambda v: v if v and str(v).strip() else 'Unknown')

    # Keep top 10 groups by count, rest becomes "Other"
    class_counts = df['Color_group'].value_counts()
    top_classes = class_counts.nlargest(10).index.tolist()
    df['Color_group'] = df['Color_group'].apply(
        lambda v: v if v in top_classes else 'Other')

    if tax_info and tax_info['subdir'] == '24_taxonomy_gtdbtk':
        rank_label = 'Class'
    else:
        rank_label = 'Genus'

    samples = sorted(df['Sample'].unique())
    # Put "Other" and "Unknown" last in the legend
    groups = sorted([g for g in df['Color_group'].unique()
                     if g not in ('Other', 'Unknown')])
    if 'Other' in df['Color_group'].values:
        groups.append('Other')
    if 'Unknown' in df['Color_group'].values:
        groups.append('Unknown')

    # Assign colors
    color_map = {}
    for i, g in enumerate(groups):
        if g == 'Other':
            color_map[g] = '#95A5A6'
        elif g == 'Unknown':
            color_map[g] = '#BDC3C7'
        else:
            color_map[g] = PALETTE[i % len(PALETTE)]

    fig = go.Figure()
    for grp in groups:
        counts = []
        for s in samples:
            counts.append(len(df[(df['Sample'] == s) & (df['Color_group'] == grp)]))
        fig.add_trace(go.Bar(
            name=grp if grp else 'Unknown', x=samples, y=counts,
            marker_color=color_map[grp]))

    fig.update_layout(
        barmode='stack', height=450,
        title=dict(text=f'Per-Sample Taxonomic Composition ({rank_label})',
                   font_color=BLUE, font_size=14),
        xaxis_tickangle=-40,
        yaxis_title='N MAGs',
        margin=dict(l=60, r=30, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.04,
                    xanchor='center', x=0.5, font_size=9),
        font=dict(family='Arial, sans-serif', size=10))
    return fig


def chart_bins_per_binner(bins_df):
    """Grouped bar chart: bins per binner per sample."""
    if bins_df.empty:
        return None
    samples = sorted(bins_df['Sample'].unique())
    binners = ['MetaBAT2', 'MaxBin2', 'SemiBin2', 'DAS Tool']
    binner_colors = [BLUE, ORANGE, GREEN, '#9B59B6']
    fig = go.Figure()
    for binner, color in zip(binners, binner_colors):
        binner_data = bins_df[bins_df['Binner'] == binner]
        counts = []
        for s in samples:
            row = binner_data[binner_data['Sample'] == s]
            counts.append(int(row['Bins'].values[0]) if not row.empty else 0)
        fig.add_trace(go.Bar(
            name=binner, x=samples, y=counts, marker_color=color,
            text=counts, textposition='outside', textfont_size=8))
    fig.update_layout(
        barmode='group', height=400,
        title=dict(text='Bins per Binner and Sample',
                   font_color=BLUE, font_size=14),
        xaxis_tickangle=-40,
        yaxis_title='N bins',
        margin=dict(l=60, r=30, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.04,
                    xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=10))
    return fig


def html_assembly_table(quast_df):
    """Table of assembly metrics from QUAST data."""
    if quast_df.empty:
        return '<p class="no-data">No assembly data (QUAST).</p>'
    rows = ''
    for _, r in quast_df.iterrows():
        sample = str(r.get('Sample', ''))
        contigs = int(r.get('# contigs (>= 0 bp)', 0))
        total_len = int(r.get('Total length (>= 0 bp)', 0))
        total_mb = total_len / 1e6
        largest = int(r.get('Largest contig', 0))
        n50 = int(r.get('N50', 0)) if 'N50' in r.index else 0
        gc = r.get('GC (%)', '-')
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{sample}">{sample}</td>'
            f'<td>{contigs:,}</td>'
            f'<td>{total_mb:,.1f}</td>'
            f'<td>{largest:,}</td>'
            f'<td>{n50:,}</td>'
            f'<td>{gc}</td>'
            f'</tr>\n')
    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Total Contigs</th><th>Total Length (Mb)</th>'
        '<th>Largest Contig</th><th>N50</th><th>GC (%)</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>')


# Colour palette for AMR-class annotation track. Stable across runs so
# the same antibiotic class always uses the same colour. Anything not
# listed falls back to a neutral grey.
AMR_CLASS_PALETTE = {
    'BETA-LACTAM'    : '#E74C3C',
    'CARBAPENEM'     : '#C0392B',
    'CEPHALOSPORIN'  : '#A93226',
    'AMINOGLYCOSIDE' : '#E67E22',
    'TETRACYCLINE'   : '#8E44AD',
    'QUINOLONE'      : '#2C5F8A',
    'FLUOROQUINOLONE': '#1F618D',
    'MACROLIDE'      : '#1ABC9C',
    'SULFONAMIDE'    : '#27AE60',
    'PHENICOL'       : '#7D6608',
    'GLYCOPEPTIDE'   : '#9C640C',
    'COLISTIN'       : '#5B2C6F',
    'TRIMETHOPRIM'   : '#117864',
    'BIOCIDE'        : '#566573',
    'STRESS'         : '#7F8C8D',
}


def _amr_class_color(cls):
    return AMR_CLASS_PALETTE.get(str(cls).upper(), '#95A5A6')


# Filters used by the gene-level AMR heatmap. Picked to match the
# review thresholds discussed in #21 — strict enough that a hit shown
# is one we would act on, lax enough that we do not lose real signal
# at sub-100% identity.
AMR_HEATMAP_MIN_IDENTITY = 90.0
AMR_HEATMAP_MIN_COVERAGE = 80.0


def chart_amr_gene_heatmap(integ_df, amr_mags_df):
    """Organism × gene heatmap with hierarchical clustering on rows and
    a coloured AMR-class annotation track on top.

    Strategy:
      • Source rows from the AMR-pathogen integration TSV (one row per
        gene call with organism, MAG, identity, coverage, location, risk).
      • Filter to high-confidence calls (identity ≥ 90%, coverage ≥ 80%).
      • Pivot to organism × gene, value = best COVERAGE per cell.
      • Cluster organisms with average linkage on euclidean distance over
        the coverage vectors (zero-imputed for missing genes); fall back
        to abundance ordering if scipy is unavailable.
      • Sort genes by AMR class then by descending number of carriers,
        so the same class clusters together visually.

    Why coverage and not identity:
      AMRFinderPlus identity tends to saturate above 95% after polished
      Nanopore assembly, so once the filter (identity ≥ 90%) is applied
      every cell falls in the same colour band — uninformative. Coverage
      keeps a real gradient: 100% = full reference gene reconstructed,
      80-90% = truncated by contig edge or partial assembly, which is
      clinically meaningful (a half-gene call is dubious).

    Honesty note: AMRFinderPlus does NOT report read depth on contig
    calls — it would be misleading to colour the heatmap by depth here.
    The KMA section is where the depth side lives.
    """
    if integ_df is None or integ_df.empty:
        return None

    df = integ_df.copy()
    df['Identity_pct'] = pd.to_numeric(df.get('Identity_pct', 0), errors='coerce').fillna(0)
    df['Coverage_pct'] = pd.to_numeric(df.get('Coverage_pct', 0), errors='coerce').fillna(0)
    df = df[(df['Identity_pct'] >= AMR_HEATMAP_MIN_IDENTITY)
            & (df['Coverage_pct'] >= AMR_HEATMAP_MIN_COVERAGE)].copy()
    if df.empty:
        return None

    df['Organism'] = df['Organism'].fillna('').astype(str).replace('', 'Unknown')
    df['MAG_label'] = df['Sample'].astype(str) + ' · ' + df['Organism']
    df['Gene'] = df['Gene'].astype(str)
    df['AMR_class'] = df['AMR_class'].astype(str).str.upper().replace('', 'OTHER')

    # Pivot: best COVERAGE per (MAG_label, Gene). Coverage is the colour
    # variable; identity stays in the hover so both signals are visible.
    pivot = df.groupby(['MAG_label', 'Gene'])['Coverage_pct'].max().unstack(
        fill_value=0)
    if pivot.empty:
        return None

    # Class lookup: most-common class for each gene (some genes appear
    # under several classes; pick the most-frequent one)
    class_lookup = (df.groupby('Gene')['AMR_class']
                      .agg(lambda x: x.value_counts().index[0])
                      .to_dict())
    # Order genes by class, then by carrier count desc, then by name
    n_carriers = (df.groupby('Gene')['MAG_label'].nunique()
                  .reindex(pivot.columns).fillna(0))
    gene_order = sorted(
        pivot.columns,
        key=lambda g: (class_lookup.get(g, 'ZZZ'),
                       -int(n_carriers.get(g, 0)),
                       g))
    pivot = pivot[gene_order]

    # Cluster rows
    row_order = list(pivot.index)
    if len(row_order) >= 3:
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            mat = pivot.values
            dists = pdist(mat, metric='euclidean')
            Z = linkage(dists, method='average')
            order = leaves_list(Z)
            row_order = [pivot.index[i] for i in order]
        except Exception:
            row_order = sorted(pivot.index,
                               key=lambda r: -float(pivot.loc[r].sum()))
    else:
        row_order = sorted(pivot.index,
                           key=lambda r: -float(pivot.loc[r].sum()))
    pivot = pivot.loc[row_order]

    # Auxiliary pivots for the hover (identity, method, location).
    # `pivot` itself now holds coverage — the colour metric.
    ident_pivot = df.groupby(['MAG_label', 'Gene'])['Identity_pct'].max().unstack(
        fill_value=0).reindex(index=pivot.index, columns=pivot.columns,
                              fill_value=0)
    method_pivot = df.groupby(['MAG_label', 'Gene'])['Method'].agg(
        lambda x: ', '.join(sorted(set(x)))).unstack(fill_value='').reindex(
            index=pivot.index, columns=pivot.columns, fill_value='')
    loc_pivot = df.groupby(['MAG_label', 'Gene'])['Location'].agg(
        lambda x: ', '.join(sorted(set(x)))).unstack(fill_value='').reindex(
            index=pivot.index, columns=pivot.columns, fill_value='')
    hover = []
    for org in pivot.index:
        row_hover = []
        for gene in pivot.columns:
            cov = pivot.loc[org, gene]
            ident = ident_pivot.loc[org, gene]
            method = method_pivot.loc[org, gene]
            loc = loc_pivot.loc[org, gene]
            cls = class_lookup.get(gene, '')
            if cov == 0:
                row_hover.append(
                    f'{org}<br>{gene} ({cls})<br>not detected at thresholds')
            else:
                row_hover.append(
                    f'{org}<br><b>{gene}</b> ({cls})<br>'
                    f'<b>coverage {cov:.1f}%</b> · identity {ident:.1f}%<br>'
                    f'method {method or "-"} · location {loc or "-"}')
        hover.append(row_hover)

    # Annotation row above the heatmap: one coloured cell per gene
    class_colors = [_amr_class_color(class_lookup.get(g, 'OTHER'))
                    for g in pivot.columns]

    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.04, 0.96], shared_xaxes=True, vertical_spacing=0.005)

    # Annotation track — single-row heatmap of class indices, with a
    # discrete colorscale built per-gene
    fig.add_trace(go.Heatmap(
        z=[[i for i in range(len(pivot.columns))]],
        x=pivot.columns.tolist(), y=['AMR class'],
        colorscale=[[i / max(1, len(pivot.columns) - 1), class_colors[i]]
                    for i in range(len(pivot.columns))],
        showscale=False,
        hovertemplate='%{x}<br>class: %{customdata}<extra></extra>',
        customdata=[[class_lookup.get(g, '') for g in pivot.columns]]),
        row=1, col=1)

    # Main heatmap. Colour = coverage %. The scale starts at the
    # filter threshold (zmin = MIN_COVERAGE) instead of 0 so the
    # full colour range is used between the meaningful 80–100% band;
    # white still encodes "not detected at thresholds" because those
    # cells are filled with 0 in the pivot and clamped below zmin.
    fig.add_trace(go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale=[
            [0.0,  '#FFFFFF'],   # below filter threshold = absent
            [0.01, '#FFFFFF'],
            [0.05, '#FDE0DD'],   # just above threshold (~80%)
            [0.30, '#FA9FB5'],
            [0.60, '#DD3497'],
            [1.0,  '#7A0177'],   # full coverage 100% = strongest
        ],
        zmin=AMR_HEATMAP_MIN_COVERAGE - 1,   # lift the floor close to the threshold
        zmax=100,
        text=hover, hovertemplate='%{text}<extra></extra>',
        colorbar=dict(title=dict(text='Coverage %', side='right'),
                      thickness=12, len=0.55, y=0.45,
                      tickvals=[AMR_HEATMAP_MIN_COVERAGE, 90, 95, 100],
                      ticktext=[f'{int(AMR_HEATMAP_MIN_COVERAGE)}%',
                                '90%', '95%', '100%'])),
        row=2, col=1)

    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(tickfont_size=9, row=1, col=1)
    fig.update_xaxes(tickangle=-50, tickfont_size=8, row=2, col=1,
                     showticklabels=True)
    fig.update_yaxes(tickfont_size=9, row=2, col=1)

    fig.update_layout(
        height=max(360, len(pivot) * 22 + 220),
        title=dict(text=('AMR profile: organisms × genes · colour = '
                         'coverage % '
                         f'(filter: ≥{AMR_HEATMAP_MIN_IDENTITY:.0f}% '
                         f'identity, ≥{AMR_HEATMAP_MIN_COVERAGE:.0f}% '
                         'coverage)'),
                   font_color=BLUE, font_size=14),
        margin=dict(l=300, r=80, t=70, b=160),
        font=dict(family='Arial, sans-serif'),
    )
    return fig


def html_amr_class_legend():
    """Static legend for the AMR-class colour annotation. Always shown
    next to the gene-level heatmap so users can decode the top track."""
    chips = []
    for cls, color in AMR_CLASS_PALETTE.items():
        chips.append(
            f'<span class="amr-class-chip" style="background:{color}">'
            f'{_esc(cls.title())}</span>')
    return ('<div class="amr-class-legend">'
            '<strong>AMR class colour key</strong> '
            + ' '.join(chips) +
            ' <span class="amr-class-chip" style="background:#95A5A6">other</span>'
            '</div>')


def chart_amr_heatmap(integ_df):
    if integ_df.empty:
        return None
    pivot = integ_df.groupby(['Organism', 'AMR_class']).size().unstack(fill_value=0)
    if pivot.empty:
        return None
    top_org = pivot.sum(axis=1).nlargest(15).index
    top_cls = pivot.sum(axis=0).nlargest(10).index
    pivot = pivot.loc[pivot.index.isin(top_org), pivot.columns.isin(top_cls)]
    if pivot.empty:
        return None
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=[[0, '#FFFFFF'], [0.01, '#F0F4F7'], [0.1, '#5B9ABF'],
                     [0.5, '#2C5F8A'], [1, '#1A3F5C']],
        text=[[str(v) if v > 0 else '' for v in row] for row in pivot.values],
        texttemplate='%{text}', textfont_size=9,
        colorbar=dict(title=dict(text='N genes', side='right'))))
    fig.update_layout(
        title=dict(text='AMR: Organisms vs Resistance Classes',
                   font_color=BLUE, font_size=14),
        height=max(350, len(pivot) * 28 + 150),
        margin=dict(l=250, r=80, t=60, b=150),
        xaxis=dict(tickangle=-45, tickfont_size=8),
        yaxis=dict(tickfont_size=9),
        font=dict(family='Arial, sans-serif'))
    return fig


def chart_risk_summary(integ_df):
    if integ_df.empty:
        return None
    order = ['CRITICAL', 'HIGH', 'MEDIUM', 'LOW']
    risk_counts = integ_df['Risk_level'].value_counts()
    # Reindex to desired order
    risk_counts = risk_counts.reindex(order).fillna(0).astype(int)
    risk_counts = risk_counts[risk_counts > 0]
    if risk_counts.empty:
        return None
    colors_map = [RISK_COLORS.get(r, GREY) for r in risk_counts.index]
    fig = go.Figure(go.Bar(
        x=risk_counts.index, y=risk_counts.values,
        marker_color=colors_map,
        text=risk_counts.values, textposition='outside'))
    fig.update_layout(
        title=dict(text='AMR Risk Level Distribution',
                   font_color=BLUE, font_size=14),
        height=300,
        margin=dict(l=60, r=30, t=60, b=60),
        yaxis_title='N genes',
        font=dict(family='Arial, sans-serif'))
    return fig


# ═══════════════════════════════════════════════════════════════
# PIPELINE DAG
# ═══════════════════════════════════════════════════════════════

def html_pipeline_dag(tax_info):
    tax_box = (f'26 {_esc(tax_info["label"])} '
               f'{_esc(tax_info["version"])} ({_esc(tax_info["db_label"])})')
    return f"""<div class="dag-container"><div class="dag-flow">
    <div class="dag-step-num">Steps 01-04</div>
    <div class="dag-row">
        <div class="dag-box qc">01 NanoPlot 1.46.2</div>
        <div class="dag-box qc">02 FastQC 0.12.1</div>
        <div class="dag-box qc">03 CHECK_SQCORE</div>
        <div class="dag-box qc">04 FastQ Screen 0.16.0</div>
    </div>
    <div class="dag-label">QC + Screening (raw reads)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 05-07</div>
    <div class="dag-row">
        <div class="dag-box trim">05 Porechop ABI 0.5.1</div>
        <div class="dag-arrow-h">&#x25B6;</div>
        <div class="dag-box trim">06 Chopper 0.12.0</div>
        <div class="dag-arrow-h">&#x25B6;</div>
        <div class="dag-box trim">07 Host Removal (minimap2 2.30)</div>
    </div>
    <div class="dag-label">Trimming + Filtering + Decontamination</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 08-09</div>
    <div class="dag-row">
        <div class="dag-box qc">08 NanoPlot 1.46.2</div>
        <div class="dag-box qc">09 FastQC 0.12.1</div>
    </div>
    <div class="dag-label">Filtered reads QC</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 10-13</div>
    <div class="dag-row">
        <div class="dag-box tax">10 Kraken2 2.17.1</div>
        <div class="dag-box tax">11 Bracken 3.1</div>
        <div class="dag-box tax">12 Kaiju 1.10.1</div>
        <div class="dag-box tax">13 Sylph 0.9.0</div>
        <div class="dag-box ann">14 KMA 1.6.8</div>
    </div>
    <div class="dag-label">EpiTax: Taxonomy + AMR reads (ResFinder)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 15-16</div>
    <div class="dag-row">
        <div class="dag-box rep">15 MultiQC 1.33</div>
        <div class="dag-box rep">16 EpiTax HTML Report</div>
    </div>
    <div class="dag-label">TAX phase reports</div>
    <div class="dag-arrow">&#x25BC; filtered reads</div>

    <div class="dag-step-num">Steps 17-19</div>
    <div class="dag-row">
        <div class="dag-box asm">17 MetaFlye 2.9.6</div>
        <div class="dag-arrow-h">&#x25B6;</div>
        <div class="dag-box asm">18 Medaka 2.2.1</div>
        <div class="dag-arrow-h">&#x25B6;</div>
        <div class="dag-box qc">19 QUAST 5.3.0</div>
    </div>
    <div class="dag-label">Assembly + Polishing + QC</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Step 20</div>
    <div class="dag-row">
        <div class="dag-box asm">20 minimap2 2.30 + samtools 1.21</div>
    </div>
    <div class="dag-label">Coverage mapping</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 21-23</div>
    <div class="dag-row">
        <div class="dag-box bin">21 MetaBAT2 2.18</div>
        <div class="dag-box bin">22 MaxBin2 2.2.7</div>
        <div class="dag-box bin">23 SemiBin2 2.2.1</div>
    </div>
    <div class="dag-label">Binning (3 algorithms)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Step 24</div>
    <div class="dag-row">
        <div class="dag-box bin">24 DAS Tool 1.1.7</div>
    </div>
    <div class="dag-label">Bin refinement</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 25-28</div>
    <div class="dag-row">
        <div class="dag-box qc">25 CheckM2 1.1.0</div>
        <div class="dag-box tax">{tax_box}</div>
        <div class="dag-box ann">27 AMRFinderPlus 4.2.7</div>
        <div class="dag-box ann">28 geNomad 1.12.0</div>
    </div>
    <div class="dag-label">QC + Taxonomy + AMR contigs + Plasmids</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Steps 29-31</div>
    <div class="dag-row">
        <div class="dag-box ann">29 ABRicate/VFDB 1.0.1</div>
        <div class="dag-box ann">30 MOB-suite 3.1.9</div>
        <div class="dag-box ann">31 IntegronFinder 2.0.5</div>
        <div class="dag-box rep">32 Bakta 1.12.0</div>
    </div>
    <div class="dag-label">Virulence + Mobility + Integrons + Annotation</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-step-num">Final step</div>
    <div class="dag-row">
        <div class="dag-box rep">AMR-Pathogen-Plasmid integration</div>
        <div class="dag-arrow-h">&#x25B6;</div>
        <div class="dag-box rep-final">EpiTaxMAG Report v2</div>
    </div>
    <div class="dag-label">Final report with risk assessment</div>
    </div></div>"""


# ═══════════════════════════════════════════════════════════════
# HTML TABLES
# ═══════════════════════════════════════════════════════════════

def html_survival_table(surv):
    if surv.empty:
        return '<p class="no-data">No QC data.</p>'
    rows = ''
    for _, r in surv.iterrows():
        pct_cls = 'good' if r['pct_reads'] >= 70 else (
            'warn' if r['pct_reads'] >= 50 else 'bad')
        gb_cls = 'highlight' if r['clean_gb'] >= 1.0 else ''
        rows += (
            f'<tr><td class="sample" data-sample="{r["Sample"]}">{r["Sample"]}</td>'
            f'<td>{r["raw_reads"]:,.0f}</td><td>{r["raw_gb"]:.2f}</td>'
            f'<td>{r["clean_reads"]:,.0f}</td><td class="{gb_cls}">{r["clean_gb"]:.2f}</td>'
            f'<td class="{pct_cls}">{r["pct_reads"]:.1f}%</td></tr>\n')
    tot = surv[['raw_reads', 'raw_gb', 'clean_reads', 'clean_gb']].sum()
    pct = (tot['clean_reads'] / tot['raw_reads'] * 100) if tot['raw_reads'] > 0 else 0
    rows += (
        f'<tr class="total-row"><td><strong>TOTAL</strong></td>'
        f'<td><strong>{tot["raw_reads"]:,.0f}</strong></td>'
        f'<td><strong>{tot["raw_gb"]:.2f}</strong></td>'
        f'<td><strong>{tot["clean_reads"]:,.0f}</strong></td>'
        f'<td><strong>{tot["clean_gb"]:.2f}</strong></td>'
        f'<td><strong>{pct:.1f}%</strong></td></tr>')
    pct_arith = float(surv['pct_reads'].mean()) if not surv.empty else 0
    note = (
        '<p class="note"><strong>The TOTAL row is read-weighted</strong> '
        f'({pct:.1f}% = sum of clean reads / sum of raw reads), which is '
        'the honest measure of run-wide data loss. The <strong>arithmetic '
        f'mean of per-sample retention</strong> is {pct_arith:.1f}% — it '
        'treats every sample equally regardless of size and tends to drop '
        'when a small low-yield sample is included. Use the weighted '
        'TOTAL for run-level reporting; the arithmetic mean is useful as '
        'a quality indicator across roughly equal samples.</p>')
    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Raw Reads</th><th>Raw Gb</th>'
        '<th>Clean Reads</th><th>Clean Gb</th><th>% Retention</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>{note}')


def _taxonomy_section_caption(tax_info):
    if tax_info['subdir'] == '24_taxonomy_gtdbtk':
        return ('Per-sample MAG distribution coloured by GTDB class. '
                f'Classifier: {tax_info["full"]}. Top 10 classes + Other.')
    if tax_info['subdir'] == '24_taxonomy_sourmash':
        return ('Per-sample MAG distribution coloured by genus from the best '
                f'Sourmash hit against {tax_info["db_label"]}. '
                f'Classifier: {tax_info["full"]}. Top 10 genera + Other.')
    # Skani is the default
    return ('Per-sample MAG distribution coloured by genus from the best '
            f'Skani hit against {tax_info["db_label"]} (ANI ≥ 80, AF ≥ 15). '
            f'Classifier: {tax_info["full"]}. Top 10 genera + Other.')


def _taxonomy_assignment_table(tax_df, tax_info):
    if tax_df.empty:
        return ''
    df = tax_df.copy()
    if tax_info['subdir'] == '24_taxonomy_gtdbtk':
        rows = ''
        for _, r in df.iterrows():
            rows += (
                f'<tr><td class="sample" data-sample="{_esc(r.get("Sample", ""))}">'
                f'{_esc(r.get("Sample", ""))}</td>'
                f'<td>{_esc(r.get("Bin", ""))}</td>'
                f'<td><em>{_esc(r.get("Organism", ""))}</em></td>'
                f'<td>{_esc(r.get("Phylum", ""))}</td>'
                f'<td>{_esc(r.get("Class", ""))}</td>'
                f'<td>{_esc(r.get("Order", ""))}</td>'
                f'<td>{_esc(r.get("Family", ""))}</td></tr>')
        return (
            '<table class="data-table sortable filterable"><thead><tr>'
            '<th>Sample</th><th>Bin</th><th>Organism</th>'
            '<th>Phylum</th><th>Class</th><th>Order</th><th>Family</th>'
            f'</tr></thead><tbody>{rows}</tbody></table>')
    # Skani / Sourmash: show ANI / AF and reference accession
    rows = ''
    for _, r in df.iterrows():
        ani = r.get('ANI', np.nan)
        af = r.get('AF', np.nan)
        ani_s = f'{ani:.2f}' if pd.notna(ani) else '-'
        af_s = f'{af:.1f}' if pd.notna(af) else '-'
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r.get("Sample", ""))}">'
            f'{_esc(r.get("Sample", ""))}</td>'
            f'<td>{_esc(r.get("Bin", ""))}</td>'
            f'<td><em>{_esc(r.get("Organism", ""))}</em></td>'
            f'<td>{_esc(r.get("Reference", ""))}</td>'
            f'<td>{ani_s}</td><td>{af_s}</td></tr>')
    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Bin</th><th>Best hit (genus species)</th>'
        '<th>Reference accession</th><th>ANI %</th><th>AF %</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>')


def html_mag_catalog(checkm2_df, tax_df, bakta_df, amr_mags_df):
    if checkm2_df.empty:
        return '<p class="no-data">No MAGs recovered.</p>'
    df = checkm2_df.copy()
    if not tax_df.empty:
        df = df.merge(
            tax_df[['Sample', 'Bin', 'Organism', 'Taxonomy']],
            on=['Sample', 'Bin'], how='left')
    else:
        df['Organism'] = ''
        df['Taxonomy'] = ''
    bakta_cols = ['tRNAs', 'rRNAs', 'CDSs', 'CRISPRs', 'ncRNAs']
    if not bakta_df.empty:
        avail = [c for c in bakta_cols if c in bakta_df.columns]
        df = df.merge(
            bakta_df[['Sample', 'Bin'] + avail].copy(),
            on=['Sample', 'Bin'], how='left')
    for c in bakta_cols:
        if c not in df.columns:
            df[c] = np.nan

    # Per-bin AMR detail: chips with gene + class colour. The integer
    # count alone was uninformative — what the reader needs to act on
    # is which genes the bin carries.
    amr_chips_per_bin = {}
    amr_count_per_bin = {}
    if amr_mags_df is not None and not amr_mags_df.empty:
        for (sample, bin_name), grp in amr_mags_df.groupby(['Sample', 'Bin']):
            amr_count_per_bin[(sample, bin_name)] = len(grp)
            chips = []
            seen = set()
            for _, r in grp.iterrows():
                gene = str(r.get('Gene', '')).strip()
                cls = str(r.get('Class', '')).strip()
                if (gene, cls) in seen or not gene:
                    continue
                seen.add((gene, cls))
                color = _amr_class_color(cls)
                chips.append(
                    f'<span class="amr-chip">'
                    f'<a href="{_card_search_url(gene)}" target="_blank" '
                    f'rel="noopener" title="CARD ARO search · {_esc(cls)}">'
                    f'<strong>{_esc(gene)}</strong></a>'
                    + (f'<span class="amr-class-chip" '
                       f'style="background:{color}">{_esc(cls.title())}</span>'
                       if cls else '')
                    + '</span>')
            amr_chips_per_bin[(sample, bin_name)] = ' '.join(chips)
    df = df.fillna('')

    def _fmt_int(v):
        try:
            return str(int(float(v))) if v != '' else '-'
        except (ValueError, TypeError):
            return '-'

    rows = ''
    for _, r in df.iterrows():
        bin_key = (r.get('Sample', ''), r.get('Bin', ''))
        amr_n = amr_count_per_bin.get(bin_key, 0)
        amr_chips = amr_chips_per_bin.get(bin_key, '')
        if amr_n > 0:
            amr_cell = (f'<span class="badge" style="background:{ORANGE}">'
                        f'{amr_n}</span>')
            amr_genes_cell = (f'<div class="mag-amr-cell">{amr_chips}</div>')
        else:
            amr_cell = '0'
            amr_genes_cell = '<span class="note">no AMR hits</span>'
        rows += (
            f'<tr><td class="sample" data-sample="{r.get("Sample", "")}">'
            f'{r.get("Sample", "")}</td>'
            f'<td>{_esc(r["Bin"])}</td>'
            f'<td><em>{_esc(r.get("Organism", ""))}</em></td>'
            f'<td>{r["Completeness"]:.1f}%</td>'
            f'<td>{r["Contamination"]:.1f}%</td>'
            f'<td>{QC_BADGE.get(r["Quality"], "")}</td>'
            f'<td>{_fmt_int(r["CDSs"])}</td>'
            f'<td>{_fmt_int(r["tRNAs"])}</td>'
            f'<td>{_fmt_int(r["rRNAs"])}</td>'
            f'<td>{_fmt_int(r["ncRNAs"])}</td>'
            f'<td>{_fmt_int(r["CRISPRs"])}</td>'
            f'<td>{amr_cell}</td>'
            f'<td>{amr_genes_cell}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Bin</th><th>Organism</th>'
        '<th>Completeness</th><th>Contamination</th><th>Quality</th>'
        '<th>CDS</th><th>tRNAs</th><th>rRNAs</th><th>ncRNAs</th>'
        '<th>CRISPR</th><th>N</th>'
        '<th style="min-width:280px">AMR genes detected (CARD link · class)</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">The <strong>N</strong> column counts AMR genes '
        'per bin; the <strong>AMR genes detected</strong> column lists '
        'them with colour-coded class so a bin\'s resistance fingerprint '
        'is readable without leaving the catalogue. Each gene name links '
        'to the CARD ARO search for that family.</p>')


def _html_amr_mags_fallback(amr_mags_df):
    """Render AMR-in-MAGs straight from AMRFinderPlus when the
    AMR-Pathogen integration table came back empty (e.g. plasmid join
    skipped or pipeline ran without that step)."""
    if amr_mags_df.empty:
        return ('<p class="no-data">AMRFinderPlus produced no AMR hits in '
                'any MAG. Reads-level hits, if any, are listed in the '
                '"AMR in reads (KMA/ResFinder)" section below.</p>')
    rows = ''
    for _, r in amr_mags_df.iterrows():
        ident = float(r.get('Identity', 0))
        cov = float(r.get('Coverage', 0))
        ident_cls = 'good' if ident >= 95 else ('warn' if ident >= 80 else 'bad')
        cov_cls = 'good' if cov >= 90 else ('warn' if cov >= 60 else 'bad')
        gene_sym = _esc(str(r.get('Gene', '')))
        accession = _esc(str(r.get('Accession', '')))
        acc_link = (f'<a href="https://www.ncbi.nlm.nih.gov/protein/{accession}" '
                    f'target="_blank" rel="noopener">{accession}</a>'
                    if accession and accession != 'NA' else '-')
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(str(r.get("Sample", "")))}">'
            f'{_esc(str(r.get("Sample", "")))}</td>'
            f'<td>{_esc(str(r.get("Bin", "")))}</td>'
            f'<td>{_esc(str(r.get("Contig", "")))}</td>'
            f'<td><strong>{gene_sym}</strong></td>'
            f'<td>{_esc(str(r.get("Class", "")))}</td>'
            f'<td>{_esc(str(r.get("Subclass", "")))}</td>'
            f'<td>{_esc(str(r.get("Method", "")))}</td>'
            f'<td class="{ident_cls}">{ident:.1f}%</td>'
            f'<td class="{cov_cls}">{cov:.1f}%</td>'
            f'<td>{acc_link}</td>'
            f'<td>{_esc(str(r.get("Closest_ref", "")))}</td></tr>')
    note = ('<p class="note">Showing raw AMRFinderPlus hits — the '
            'plasmid/chromosome integration table returned no rows for this run, '
            'so confidence scoring is not applied here.</p>')
    return note + (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>MAG</th><th>Contig</th><th>Gene</th>'
        '<th>Class</th><th>Subclass</th><th>Method</th>'
        '<th>Identity</th><th>Coverage</th><th>Accession</th>'
        '<th>Closest reference</th></tr></thead>'
        f'<tbody>{rows}</tbody></table>')


def html_amr_detail(integ_df, amr_mags_df):
    if integ_df.empty:
        return _html_amr_mags_fallback(amr_mags_df)
    df = integ_df.copy()

    # Build lookup from amr_mags for extra columns
    amr_lookup = {}
    if not amr_mags_df.empty:
        for _, ar in amr_mags_df.iterrows():
            key = (str(ar.get('Sample', '')), str(ar.get('Bin', '')),
                   str(ar.get('Gene', '')), str(ar.get('Contig', '')))
            amr_lookup[key] = {
                'Accession': _safe_str(ar.get('Accession', '')),
                'Closest_ref': _safe_str(ar.get('Closest_ref', '')),
                'HMM_desc': _safe_str(ar.get('HMM_desc', '')),
            }

    rows = ''
    for _, r in df.iterrows():
        conf = confidence_score(
            r.get('Identity_pct', 0), r.get('Coverage_pct', 0),
            r.get('Method', ''),
            r.get('Location', '') == 'PLASMID',
            r.get('MAG_quality', ''),
            r.get('Plasmid_score', 0))
        conf_cls = 'good' if conf >= 70 else ('warn' if conf >= 40 else 'bad')
        loc = r.get('Location', '')
        if loc == 'PLASMID':
            loc_badge = f'<span class="badge" style="background:{RED}">PLASMID</span>'
        else:
            loc_badge = f'<span class="badge" style="background:{LBLUE}">CHROMOSOME</span>'
        s_name = str(r.get('Sample', ''))

        gene_sym = str(r.get('Gene', ''))
        mag_name = str(r.get('MAG', r.get('Bin', '')))
        contig_name = str(r.get('Contig', ''))
        lookup_key = (s_name, mag_name, gene_sym, contig_name)
        extra = amr_lookup.get(lookup_key, {})
        accession = extra.get('Accession', '')
        closest_ref = extra.get('Closest_ref', '')
        if not accession and not closest_ref:
            for k, v in amr_lookup.items():
                if k[0] == s_name and k[2] == gene_sym:
                    accession = v.get('Accession', '')
                    closest_ref = v.get('Closest_ref', '')
                    break

        description = closest_ref if closest_ref else str(
            r.get('Gene_description', ''))

        links_parts = []
        if accession:
            ncbi_url = f'https://www.ncbi.nlm.nih.gov/protein/{accession}'
            links_parts.append(
                f'<a href="{ncbi_url}" target="_blank" rel="noopener" '
                f'title="NCBI Protein record for the closest reference">'
                f'{_esc(accession)}</a>')
        if gene_sym:
            card_url = _card_search_url(gene_sym)
            links_parts.append(
                f'<a href="{card_url}" target="_blank" rel="noopener" '
                f'title="CARD ARO search · {_esc(gene_sym)}">CARD</a>')
        links_html = ' | '.join(links_parts) if links_parts else '-'

        method_val = str(r.get('Method', ''))

        risk_rule = str(r.get('Risk_rule', '') or '')
        risk_reason = str(r.get('Risk_reason', '') or '')
        risk_cell = (
            f'{RISK_BADGE.get(r.get("Risk_level", ""), "")} '
            f'<small class="risk-rule" title="{_esc(risk_reason)}">'
            f'{_esc(risk_rule)}</small>'
            if risk_rule else RISK_BADGE.get(r.get('Risk_level', ''), ''))

        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{s_name}">{s_name}</td>'
            f'<td><em>{_esc(r.get("Organism", ""))}</em></td>'
            f'<td><strong>{_esc(gene_sym)}</strong></td>'
            f'<td class="desc-cell">{_esc(description)}</td>'
            f'<td>{_esc(r.get("AMR_class", ""))}</td>'
            f'<td>{loc_badge}</td>'
            f'<td>{r.get("Identity_pct", 0):.1f}%</td>'
            f'<td>{r.get("Coverage_pct", 0):.1f}%</td>'
            f'<td>{_esc(method_val)}</td>'
            f'<td class="{conf_cls}"><strong>{conf}</strong></td>'
            f'<td>{risk_cell}</td>'
            f'<td class="links-cell">{links_html}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Organism</th><th>Gene</th><th>Description</th>'
        '<th>AMR Class</th><th>Location</th>'
        '<th>% Identity</th><th>% Coverage</th><th>Method</th>'
        '<th>Confidence</th><th>Risk (rule)</th><th>Links</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Confidence: composite score (0-100) based on identity, '
        'coverage, detection method, KMA-MAG concordance and MAG quality. '
        'Risk badge is followed by the matched rule ID (R1–R5); hover for the '
        'rule trigger or open the “Risk model” section. '
        'Links: <strong>NCBI Protein</strong> (closest reference accession) and '
        '<strong>CARD</strong> (Comprehensive Antibiotic Resistance Database — '
        'ARO search by gene name, returns variants and ontology). '
        'NCBI Pathogens RefGene was tried first but its search was '
        'unreliable for short gene families like <code>fos</code> or '
        '<code>tet</code> — CARD handles these family-level queries.</p>')


KMA_DEPTH_TRUSTED = 10.0   # ≥10× is the threshold of credible read evidence


def _render_kma_rows(df):
    rows = ''
    for _, r in df.iterrows():
        ident = float(r['KMA_Identity'])
        cov = float(r['KMA_Coverage'])
        depth = float(r['KMA_Depth'])
        ident_cls = 'good' if ident >= 95 else ('warn' if ident >= 80 else 'bad')
        cov_cls = 'good' if cov >= 90 else ('warn' if cov >= 60 else 'bad')
        depth_cls = ('good' if depth >= 30 else
                     ('warn' if depth >= KMA_DEPTH_TRUSTED else 'bad'))
        gene = str(r['Gene_root'])
        gene_link = (
            f'<a href="{_card_search_url(gene)}" '
            f'target="_blank" rel="noopener" title="CARD ARO search">'
            f'<strong>{_esc(gene)}</strong></a>')
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(str(r["Sample"]))}">'
            f'{_esc(str(r["Sample"]))}</td>'
            f'<td>{gene_link}</td>'
            f'<td><code>{_esc(str(r["Template"]))}</code></td>'
            f'<td class="{ident_cls}">{ident:.2f}%</td>'
            f'<td class="{cov_cls}">{cov:.2f}%</td>'
            f'<td class="{depth_cls}">{depth:.1f}x</td></tr>\n')
    return rows


def html_kma_reads_table(kma_df):
    """Render KMA/ResFinder hits split into trusted (≥10× depth) and
    low-depth artefact-prone tiers. The 10× threshold is conservative for
    Nanopore — at 1× the alignment is barely better than chance and any
    AMR call should be treated as non-evidence."""
    if kma_df.empty:
        return ('<p class="no-data">KMA found no AMR hits in the filtered '
                'reads (ResFinder DB).</p>')
    df = kma_df.copy()
    df['KMA_Identity'] = pd.to_numeric(df['KMA_Identity'], errors='coerce').fillna(0)
    df['KMA_Coverage'] = pd.to_numeric(df['KMA_Coverage'], errors='coerce').fillna(0)
    df['KMA_Depth'] = pd.to_numeric(df['KMA_Depth'], errors='coerce').fillna(0)
    df = df.sort_values(['Sample', 'KMA_Depth', 'Gene_root'],
                        ascending=[True, False, True])
    trusted = df[df['KMA_Depth'] >= KMA_DEPTH_TRUSTED]
    low = df[df['KMA_Depth'] < KMA_DEPTH_TRUSTED]

    def _table(sub, caption):
        if sub.empty:
            return f'<p class="no-data">{caption}: no rows.</p>'
        return (
            f'<h4 style="margin-top:1em">{caption} ({len(sub)} rows)</h4>'
            '<table class="data-table sortable filterable"><thead><tr>'
            '<th>Sample</th><th>Gene</th><th>Template (ResFinder)</th>'
            '<th>Identity</th><th>Coverage</th><th>Depth</th>'
            f'</tr></thead><tbody>{_render_kma_rows(sub)}</tbody></table>')

    return (
        _table(trusted,
               f'Trusted hits (depth ≥ {KMA_DEPTH_TRUSTED:.0f}× — credible read evidence)') +
        '<details class="prov-panel" style="margin-top:0.6rem">'
        '<summary><span class="prov-icon">+</span>'
        '<span>Low-depth tier (drill-down — likely artefacts on Nanopore)</span>'
        '</summary>'
        f'{_table(low, f"Low-depth hits (depth < {KMA_DEPTH_TRUSTED:.0f}×)")}'
        '</details>'
        '<p class="note">KMA aligns filtered reads against the ResFinder '
        'AMR template DB. On Nanopore, depths under '
        f'{KMA_DEPTH_TRUSTED:.0f}× are short-template artefacts in the vast '
        'majority of cases — high <em>coverage</em> with low <em>depth</em> '
        'is the hallmark. A trusted hit needs depth, identity and coverage '
        'all to align. Gene names link to CARD (ARO search).</p>')


_FQS_CATEGORY_BADGES = {
    'target'         : ('Target', GREEN),
    'host'           : ('Host', RED),
    'control'        : ('Control', GREY),
    'parasite'       : ('Parasite', '#9B59B6'),
    'virus'          : ('Virus', '#16A085'),
    'other_reference': ('Other ref.', LBLUE),
}


def html_fastqscreen_summary(fqs_detail_df, target_organisms):
    """Executive view: one row per sample.
    Columns: reads, %Hit_no_genomes (from file footer), best target hit,
    strongest host signal, strongest control signal, verdict.

    The "verdict" is an inference from the report; same provenance rule
    as the detail table — flagged in the footnote."""
    if fqs_detail_df is None or fqs_detail_df.empty:
        return ('<p class="no-data">No FastQ Screen detail tables found '
                '(02_fastqscreen_raw/*_screen.txt missing).</p>')
    target_set = [g.strip() for g in (target_organisms or '').split(',') if g.strip()]
    no_hit_lookup = fqs_detail_df.attrs.get('no_hit_pct', {})

    df = fqs_detail_df.copy()
    df['Category'] = df['Genome'].apply(
        lambda g: classify_panel_genome(g, target_set)[0])

    # Build per-sample rows, then sort by triage priority before rendering:
    # Contamination first, then Watch, then Clean. Within each verdict
    # bucket, the sample with the strongest contamination signal goes on
    # top so the eye lands on the most-actionable rows first.
    sample_rows = []
    for sample in sorted(df['Sample'].unique()):
        sub = df[df['Sample'] == sample]
        if sub.empty:
            continue
        reads_processed = int(sub['Reads_processed'].iloc[0])
        no_hit_pct = no_hit_lookup.get(sample)

        def best_in(category, metric):
            cs = sub[sub['Category'] == category]
            if cs.empty or cs[metric].max() == 0:
                return None, 0.0
            row = cs.loc[cs[metric].idxmax()]
            return row['Genome'], float(row[metric])

        target_g, target_v = best_in('target', 'Unique_pct')
        host_g, host_v     = best_in('host', 'Total_hit_pct')
        ctrl_g, ctrl_v     = best_in('control', 'Total_hit_pct')

        if host_v >= 1.0 or ctrl_v >= 5.0:
            verdict_color, verdict_label, verdict_rank = RED, 'Contamination', 0
        elif host_v >= 0.1 or ctrl_v >= 1.0:
            verdict_color, verdict_label, verdict_rank = ORANGE, 'Watch', 1
        else:
            verdict_color, verdict_label, verdict_rank = GREEN, 'Clean', 2

        sample_rows.append({
            'sample': sample, 'reads': reads_processed,
            'no_hit_pct': no_hit_pct,
            'target_g': target_g, 'target_v': target_v,
            'host_g': host_g, 'host_v': host_v,
            'ctrl_g': ctrl_g, 'ctrl_v': ctrl_v,
            'verdict_color': verdict_color,
            'verdict_label': verdict_label,
            'verdict_rank': verdict_rank,
            # Triage priority: strongest of host vs control
            'priority_signal': max(host_v, ctrl_v),
        })

    sample_rows.sort(key=lambda r: (r['verdict_rank'], -r['priority_signal'],
                                    r['sample']))

    summary_rows = ''
    for r in sample_rows:
        target_cell = (
            f'<strong>{_esc(r["target_g"])}</strong> · {r["target_v"]:.2f}%'
            if r['target_g'] else '<span class="note">none above 0%</span>')
        host_cell = (
            f'{_esc(r["host_g"])} · {r["host_v"]:.2f}%' if r['host_g']
            else '<span class="note">none</span>')
        ctrl_cell = (
            f'{_esc(r["ctrl_g"])} · {r["ctrl_v"]:.2f}%' if r['ctrl_g']
            else '<span class="note">none</span>')
        no_hit_cell = (f'{r["no_hit_pct"]:.2f}%' if r['no_hit_pct'] is not None
                       else '<span class="note">—</span>')

        summary_rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["sample"])}">'
            f'<strong>{_esc(r["sample"])}</strong></td>'
            f'<td>{r["reads"]:,}</td>'
            f'<td>{no_hit_cell}</td>'
            f'<td>{target_cell}</td>'
            f'<td>{host_cell}</td>'
            f'<td>{ctrl_cell}</td>'
            f'<td><span class="badge" style="background:{r["verdict_color"]}">'
            f'{r["verdict_label"]}</span></td></tr>')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Reads (subset)</th>'
        '<th>% Reads not hitting any panel genome</th>'
        '<th>Best target (Unique%)</th>'
        '<th>Strongest host signal</th>'
        '<th>Strongest control / vector signal</th>'
        '<th>Verdict</th>'
        f'</tr></thead><tbody>{summary_rows}</tbody></table>'
        '<p class="note"><strong>Rows sorted by triage priority</strong>: '
        'Contamination → Watch → Clean, then by strongest host/control '
        'signal within each bucket — not alphabetically. '
        '<strong>From the file</strong>: Reads, <code>%Hit_no_genomes</code> '
        'footer, and the per-genome percentages used to pick the “strongest” '
        'signal. <strong>Inferred by this report</strong>: the per-genome '
        'category (target / host / control) and the Verdict badge '
        '(thresholds: host ≥ 1% or controls ≥ 5% → Contamination; '
        'host ≥ 0.1% or controls ≥ 1% → Watch; otherwise Clean). '
        'Drill-down per sample × panel genome below.</p>')


def html_fastqscreen_detail(fqs_detail_df, target_organisms):
    """Per-sample × genome FastQ Screen breakdown.

    Five percentage columns come straight from `*_screen.txt` (One-hit-one,
    Multi-hit-one, One-hit-multi, Multi-hit-multi, Unmapped); the
    Category column is an inference of this report based on
    `params.target_organisms` and a panel-name lookup, NOT something
    FastQ Screen produces. The caption above the table makes that split
    explicit.
    """
    if fqs_detail_df.empty:
        return ('<p class="no-data">No FastQ Screen detail tables found '
                '(02_fastqscreen_raw/*_screen.txt missing).</p>')
    df = fqs_detail_df.copy()
    target_set = [g.strip() for g in (target_organisms or '').split(',') if g.strip()]
    df['Category'] = df['Genome'].apply(
        lambda g: classify_panel_genome(g, target_set)[0])
    # Order: target → host → control → parasite/virus → other
    cat_order = ['target', 'host', 'control', 'parasite', 'virus', 'other_reference']
    df['_cat_rank'] = df['Category'].apply(
        lambda c: cat_order.index(c) if c in cat_order else len(cat_order))
    df = df.sort_values(['Sample', '_cat_rank', 'Total_hit_pct'],
                        ascending=[True, True, False])

    rows = ''
    for _, r in df.iterrows():
        cat_label, cat_color = _FQS_CATEGORY_BADGES.get(
            r['Category'], ('?', GREY))
        cat_badge = (f'<span class="badge" style="background:{cat_color}">'
                     f'{cat_label}</span>')

        # Visual highlight of "is this read alignment trustworthy?"
        unique = r['Unique_pct']
        cross = r['CrossMapping_pct']
        repeat = r['RepeatLike_pct']
        unique_cls = ('good' if unique >= 1.0 else
                      ('warn' if unique >= 0.1 else 'note'))
        cross_cls = ('warn' if cross >= 5 else 'note')

        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["Sample"])}">'
            f'{_esc(r["Sample"])}</td>'
            f'<td><strong>{_esc(r["Genome"])}</strong></td>'
            f'<td>{cat_badge}</td>'
            f'<td>{r["Total_hit_pct"]:.2f}%</td>'
            f'<td class="{unique_cls}"><strong>{unique:.2f}%</strong></td>'
            f'<td>{r["MultiInternal_pct"]:.2f}%</td>'
            f'<td class="{cross_cls}">{cross:.2f}%</td>'
            f'<td>{repeat:.2f}%</td>'
            f'<td>{r["Unmapped_pct"]:.2f}%</td></tr>')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Panel genome</th><th>Category</th>'
        '<th>Total hit %</th><th>Unique %</th><th>Multi (1-genome) %</th>'
        '<th>Cross-mapping %</th><th>Repeat-like %</th><th>Unmapped %</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note"><strong>From the file</strong> (FastQ Screen 0.16.0): '
        'the five percentage columns and the Reads_processed field. '
        '<strong>Inferred by this report</strong>: the Category badge '
        '(target / host / control / parasite / virus / other reference), '
        'derived from <code>params.target_organisms</code> plus a small '
        'known-panel lookup (<code>_PANEL_GENUS_LOOKUP</code>). '
        'A high <em>Unique</em> % is the trustworthy signal — multi-hit-one '
        'reflects internal homology of the reference, '
        'cross-mapping is shared content with another panel reference, and '
        'repeat-like is reads with multiple hits across multiple references '
        '(noisy signal, often low-complexity or ambiguous reads).</p>'
        '<p class="note"><strong>Why high cross-mapping / repeat-like '
        'often appear (&gt;30%)</strong>: the curated panel includes '
        'several closely related references (multiple Acinetobacter, '
        'multiple Bacteroides, several Enterobacteriaceae, two OXA-bearing '
        'Acinetobacter spp.). When a read comes from a real organism in '
        'the sample, it tends to map to <em>several</em> panel members at '
        'once because they share &gt;90% nucleotide identity over large '
        'regions. That is not contamination — it is the natural cost of a '
        'phylogenetically dense panel for surveillance. The honest '
        'signals to act on are: <em>Unique %</em> (clean target detection), '
        '<em>Host %</em> at the host row (real contamination), and '
        '<em>Control %</em> at PhiX/vector rows (lab contamination). '
        'Cross-mapping concentrated within a phylogenetically tight subset '
        'of the panel (e.g. all Acinetobacter rows lit) is a positive '
        'signal, not a problem.</p>')


def html_sylph_summary(comparison_df):
    """Executive view: one row per sample with shared/discordant counts,
    top-3 disagreements (genus + which side detected it + abundance), and
    a single verdict badge.

    Rows sorted by total disagreements descending — samples where the
    two classifiers part ways most strongly land on top. Confirmed
    samples (zero disagreements) drop to the bottom."""
    if comparison_df is None or comparison_df.empty:
        return ''
    sample_rows = []
    for sample in sorted(comparison_df['Sample'].unique()):
        sub = comparison_df[comparison_df['Sample'] == sample]
        n_both = int((sub['Detected_by'] == 'Both').sum())
        n_b_only = int((sub['Detected_by'] == 'Bracken only').sum())
        n_s_only = int((sub['Detected_by'] == 'Sylph only').sum())
        n_total = len(sub)
        if n_total == 0:
            continue
        # Verdict + colour
        if n_s_only == 0 and n_b_only == 0:
            verdict_label = 'Confirmed'
            verdict_color = GREEN
        elif n_s_only > n_b_only:
            verdict_label = 'Sylph adds'
            verdict_color = LBLUE
        elif n_b_only > n_s_only:
            verdict_label = 'Sylph misses'
            verdict_color = ORANGE
        else:
            verdict_label = 'Symmetric'
            verdict_color = GREY

        disagree = sub[sub['Detected_by'] != 'Both'].copy()
        if not disagree.empty:
            disagree['_max_pct'] = disagree[['Bracken_pct', 'Sylph_pct']].max(axis=1)
            top3 = disagree.nlargest(3, '_max_pct')
            top3_str = ' · '.join(
                f'<em>{_esc(r["Genus"])}</em> '
                f'<small>({"B" if r["Detected_by"] == "Bracken only" else "S"} '
                f'{max(r["Bracken_pct"], r["Sylph_pct"]):.1f}%)</small>'
                for _, r in top3.iterrows())
            max_disagreement_pct = float(disagree['_max_pct'].max())
        else:
            top3_str = '<span class="note">no disagreements</span>'
            max_disagreement_pct = 0.0

        sample_rows.append({
            'sample': sample,
            'n_both': n_both, 'n_b_only': n_b_only, 'n_s_only': n_s_only,
            'top3_str': top3_str,
            'verdict_color': verdict_color, 'verdict_label': verdict_label,
            'priority': (n_b_only + n_s_only, max_disagreement_pct),
        })

    # Triage priority: most discordant samples first, ties broken by the
    # largest abundance gap and finally by sample name.
    sample_rows.sort(
        key=lambda r: (-r['priority'][0], -r['priority'][1], r['sample']))

    rows = ''
    for r in sample_rows:
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["sample"])}">'
            f'<strong>{_esc(r["sample"])}</strong></td>'
            f'<td>{r["n_both"]}</td>'
            f'<td>{r["n_b_only"]}</td>'
            f'<td>{r["n_s_only"]}</td>'
            f'<td>{r["top3_str"]}</td>'
            f'<td><span class="badge" style="background:{r["verdict_color"]}">'
            f'{r["verdict_label"]}</span></td></tr>')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Shared genera</th>'
        '<th>Bracken-only</th><th>Sylph-only</th>'
        '<th>Top-3 disagreements (B/S = which side, abundance)</th>'
        '<th>Verdict</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note"><strong>Rows sorted by total disagreements</strong> '
        '(highest first), ties broken by the largest single-genus '
        'abundance gap — Confirmed samples sink to the bottom. '
        '<strong>From the data</strong>: shared / Bracken-only / Sylph-only '
        'counts and the per-genus abundances. '
        '<strong>Inferred by this report</strong>: the Verdict badge '
        '(<strong>Confirmed</strong> = top-N agree fully; '
        '<strong>Sylph adds</strong> = Sylph detects more genera Bracken missed; '
        '<strong>Sylph misses</strong> = Bracken detects more genera Sylph missed; '
        '<strong>Symmetric</strong> = same number of disagreements on each side). '
        'Open the per-genus table below for the full audit.</p>')


def html_sylph_concordance(comparison_df):
    """Render the Bracken-vs-Sylph genus concordance table with an
    explicit honest summary at the top: does Sylph add signal here?"""
    if comparison_df.empty:
        return ('<p class="no-data">No Bracken or Sylph genus tables found '
                'for this run, so no honest concordance can be reported.</p>')

    summary_blocks = []
    for sample in sorted(comparison_df['Sample'].unique()):
        sub = comparison_df[comparison_df['Sample'] == sample]
        n_both = int((sub['Detected_by'] == 'Both').sum())
        n_b_only = int((sub['Detected_by'] == 'Bracken only').sum())
        n_s_only = int((sub['Detected_by'] == 'Sylph only').sum())
        n_total = len(sub)
        if n_total == 0:
            continue
        # Honest sentence
        if n_s_only == 0 and n_b_only == 0:
            verdict = ('<span class="good"><strong>Sylph confirms Bracken '
                       'top-N — no additional or missing genera.</strong></span>')
        elif n_s_only > n_b_only:
            verdict = (f'<strong>Sylph adds {n_s_only} genera</strong> not '
                       f'in Bracken top-N (and is missing {n_b_only}). '
                       'Plausibly real signal — check the Sylph-only column '
                       'before relying on it.')
        elif n_b_only > n_s_only:
            verdict = (f'<strong>Sylph misses {n_b_only} Bracken genera</strong> '
                       f'(and adds {n_s_only}). Common when the Sylph DB '
                       'sketch is sparser than the Kraken2 DB at low coverage.')
        else:
            verdict = (f'Symmetric difference: {n_b_only} Bracken-only · '
                       f'{n_s_only} Sylph-only · {n_both} shared. '
                       'No clear winner — keep both as orthogonal evidence.')
        summary_blocks.append(
            f'<p><code>{_esc(sample)}</code>: {verdict} '
            f'({n_both} shared / {n_b_only} Bracken-only / '
            f'{n_s_only} Sylph-only).</p>')

    rows = ''
    df = comparison_df.copy().sort_values(
        ['Sample', 'Detected_by', 'Bracken_pct', 'Sylph_pct'],
        ascending=[True, True, False, False])
    for _, r in df.iterrows():
        det = r['Detected_by']
        det_color = (GREEN if det == 'Both' else
                     ORANGE if det == 'Bracken only' else LBLUE)
        b_pct = r['Bracken_pct']
        s_pct = r['Sylph_pct']
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["Sample"])}">'
            f'{_esc(r["Sample"])}</td>'
            f'<td><strong>{_esc(r["Genus"])}</strong></td>'
            f'<td>{b_pct:.3f}%</td>'
            f'<td>{s_pct:.3f}%</td>'
            f'<td><span class="badge" style="background:{det_color}">'
            f'{det}</span></td></tr>')

    return (
        ''.join(summary_blocks) +
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Genus</th>'
        '<th>Bracken %</th><th>Sylph %</th><th>Concordance</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note">Top-15 most abundant genera per tool are compared. '
        'Bracken values come from <code>08_tax_bracken/*.bracken.G.txt</code> '
        '(<code>fraction_total_reads</code>); Sylph values come from the '
        'Taxonomic_abundance column of <code>10_tax_sylph/sylph_profile_all.tsv</code>, '
        'with the genus parsed best-effort from the reference Contig_name. '
        'A genus shown in only one column is not necessarily wrong — it is '
        'an evidence asymmetry that warrants a second look.</p>')


def html_mobsuite_detail(mobsuite_df, integ_df, amr_mags_df, tax_df):
    """Per-plasmid MOB-suite typing with cross-link to AMR carriage and
    MAG taxonomy. Honest when MOB-suite did not call any plasmid."""
    if mobsuite_df is None or mobsuite_df.empty:
        return ('<p class="no-data"><strong>MOB-suite did not call any '
                'circular plasmid for this run.</strong> '
                'mob_recon needs marker-rich, well-assembled circular contigs; '
                'metagenomic Nanopore assemblies often fragment plasmids '
                'too much for MOB-suite to close them. The geNomad table '
                'above (which scores any contig that looks plasmidic, even '
                'fragmented) remains the primary plasmid evidence.</p>')

    # AMR genes-per-contig from amr_mags + per-plasmid integration
    amr_by_contig = {}
    if amr_mags_df is not None and not amr_mags_df.empty:
        for _, r in amr_mags_df.iterrows():
            key = (str(r.get('Sample', '')), str(r.get('Contig', '')))
            amr_by_contig.setdefault(key, []).append(
                f'{r.get("Gene", "")} ({r.get("Class", "")})')
    plasmid_amr_by_sample = {}
    if integ_df is not None and not integ_df.empty:
        for _, r in integ_df[integ_df['Location'] == 'PLASMID'].iterrows():
            key = (str(r.get('Sample', '')), str(r.get('Contig', '')))
            plasmid_amr_by_sample.setdefault(key, []).append(
                f'{r.get("Gene", "")} ({r.get("AMR_class", "")})')

    rows = ''
    for _, r in mobsuite_df.iterrows():
        amr_list = sorted(set(
            amr_by_contig.get((r['Sample'], r['sample_id']), []) +
            plasmid_amr_by_sample.get((r['Sample'], r['sample_id']), [])))
        amr_html = (', '.join(amr_list) if amr_list
                    else '<span class="note">none on this plasmid</span>')
        mobility_color = (RED if r['predicted_mobility'] == 'conjugative'
                          else (ORANGE if r['predicted_mobility'] == 'mobilizable'
                                else GREY))
        rows += (
            f'<tr><td class="sample" data-sample="{_esc(r["Sample"])}">'
            f'{_esc(r["Sample"])}</td>'
            f'<td><code>{_esc(r["sample_id"])}</code></td>'
            f'<td>{int(r["size"]):,} bp · GC {r["gc"]:.1f}%</td>'
            f'<td>{_esc(r["rep_type"])}</td>'
            f'<td>{_esc(r["relaxase_type"])}</td>'
            f'<td><span class="badge" style="background:{mobility_color}">'
            f'{_esc(r["predicted_mobility"]) or "-"}</span></td>'
            f'<td>{_esc(r["primary_cluster_id"])}</td>'
            f'<td class="desc-cell">{_esc(r["mash_neighbor_identification"])}</td>'
            f'<td>{amr_html}</td></tr>')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Plasmid (mob_recon contig)</th>'
        '<th>Size</th><th>Replicon (rep)</th>'
        '<th>Relaxase (MOB)</th><th>Mobility</th>'
        '<th>Cluster</th><th>Closest reference (Mash)</th>'
        '<th>AMR genes on this plasmid</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note"><strong>Mobility</strong>: <em>conjugative</em> = '
        'self-transferable plasmid (highest dissemination risk); '
        '<em>mobilizable</em> = needs another plasmid in trans; '
        '<em>non-mobilizable</em> = no MOB markers found. The AMR-genes '
        'column is a cross-link of this report into '
        '<code>25_amr_mags/</code> and the AMR-pathogen integration TSV — '
        'shown here to make horizontal-transfer risk evident at a glance.</p>')


def html_amr_on_plasmids(amr_plas_df):
    """Render the AMRFinderPlus hits on geNomad plasmid sequences with
    a per-contig summary first (one row per contig with all its genes
    chipped together) and the per-hit detail in a collapsed drill-down.
    Multiple hits on the same contig are common — this layout makes
    co-localisation obvious without forcing a long table."""
    if amr_plas_df is None or amr_plas_df.empty:
        return ('<p class="no-data">AMRFinderPlus found no AMR genes on '
                'the geNomad-extracted plasmid contigs of this run. '
                'The geNomad <code>amr_genes</code> column above is the '
                'only AMR-on-plasmid evidence available.</p>')
    df = amr_plas_df.copy()
    df['Identity'] = pd.to_numeric(df.get('Identity', 0), errors='coerce').fillna(0)
    df['Coverage'] = pd.to_numeric(df.get('Coverage', 0), errors='coerce').fillna(0)

    # ── Per-contig summary ──
    summary_rows = ''
    for (sample, contig), grp in df.groupby(['Sample', 'Contig']):
        n_hits = len(grp)
        chips = []
        for _, hit in grp.iterrows():
            gene = str(hit.get('Gene', '')).strip()
            cls = str(hit.get('Class', '')).strip()
            ident = float(hit.get('Identity', 0))
            color = _amr_class_color(cls)
            chip = (
                f'<span class="amr-chip">'
                f'<a href="{_card_search_url(gene)}" target="_blank" rel="noopener" '
                f'title="CARD ARO search · identity {ident:.1f}%">'
                f'<strong>{_esc(gene)}</strong></a>'
                f'<span class="amr-class-chip" style="background:{color}">'
                f'{_esc(cls.title())}</span></span>')
            chips.append(chip)
        summary_rows += (
            f'<tr><td class="sample" data-sample="{_esc(str(sample))}">'
            f'{_esc(str(sample))}</td>'
            f'<td><code>{_esc(str(contig))}</code></td>'
            f'<td>{n_hits}</td>'
            f'<td class="desc-cell">{" ".join(chips)}</td></tr>')

    summary_table = (
        '<h4>Per-plasmid-contig summary · co-localised AMR genes</h4>'
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Plasmid contig</th>'
        '<th>N AMR hits</th>'
        '<th>Genes detected (CARD link · class)</th>'
        f'</tr></thead><tbody>{summary_rows}</tbody></table>'
        '<p class="note">Co-localised genes on the same plasmid contig '
        'travel together — that is the metric of horizontal-transfer '
        'risk this view surfaces. A contig carrying e.g. a β-lactamase '
        '+ a fluoroquinolone resistance gene is a higher dissemination '
        'concern than two contigs each carrying one.</p>')

    # ── Per-hit detail (drill-down) ──
    df = df.sort_values(['Sample', 'Contig', 'Identity'],
                        ascending=[True, True, False])
    detail_rows = ''
    for _, r in df.iterrows():
        ident = float(r['Identity'])
        cov = float(r['Coverage'])
        ident_cls = 'good' if ident >= 95 else ('warn' if ident >= 80 else 'bad')
        cov_cls = 'good' if cov >= 90 else ('warn' if cov >= 60 else 'bad')
        gene = str(r.get('Gene', '')).strip()
        cls = str(r.get('Class', '')).strip()
        sub = str(r.get('Subclass', '')).strip()
        gene_link = (f'<a href="{_card_search_url(gene)}" target="_blank" '
                     f'rel="noopener" title="CARD ARO search">'
                     f'<strong>{_esc(gene)}</strong></a>')
        cls_chip = (f'<span class="amr-class-chip" '
                    f'style="background:{_amr_class_color(cls)}">'
                    f'{_esc(cls.title())}</span>') if cls else ''
        detail_rows += (
            f'<tr><td class="sample" data-sample="{_esc(str(r["Sample"]))}">'
            f'{_esc(str(r["Sample"]))}</td>'
            f'<td><code>{_esc(str(r["Contig"]))}</code></td>'
            f'<td>{gene_link}</td>'
            f'<td>{cls_chip}</td>'
            f'<td><small>{_esc(sub.replace("/", " / "))}</small></td>'
            f'<td class="{ident_cls}">{ident:.1f}%</td>'
            f'<td class="{cov_cls}">{cov:.1f}%</td>'
            f'<td>{_esc(str(r.get("Method", "")))}</td>'
            f'<td class="desc-cell">{_esc(str(r.get("Name", "")))}</td>'
            f'<td>{_esc(str(r.get("Accession", "")))}</td></tr>\n')
    detail_table = (
        '<details class="prov-panel" style="margin-top:0.6rem">'
        '<summary><span class="prov-icon">+</span>'
        f'<span>Drill-down: per-hit detail ({len(df)} AMR calls)</span>'
        '</summary>'
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Plasmid contig</th><th>Gene</th>'
        '<th>Class</th><th>Subclass</th>'
        '<th>Identity</th><th>Coverage</th><th>Method</th>'
        '<th>Closest reference description</th><th>Accession</th>'
        f'</tr></thead><tbody>{detail_rows}</tbody></table>'
        '</details>')

    return (summary_table + detail_table +
            '<p class="note">Independent AMRFinderPlus pass on the FASTA '
            'of plasmid contigs predicted by geNomad — orthogonal to the '
            '<code>amr_genes</code> column shown above.</p>')


def html_concordance_table(conc_df):
    if conc_df.empty:
        return '<p class="no-data">No KMA-MAG concordance data.</p>'
    rows = ''
    for _, r in conc_df.iterrows():
        badge = CONCORDANCE_BADGE.get(r.get('Source', ''), '')
        kma_d = (f"{r.get('KMA_Depth', 0):.1f}x"
                 if pd.notna(r.get('KMA_Depth')) and r.get('KMA_Depth', 0) > 0
                 else '-')
        kma_i = (f"{r.get('KMA_Identity', 0):.1f}%"
                 if pd.notna(r.get('KMA_Identity')) and r.get('KMA_Identity', 0) > 0
                 else '-')
        mag_i = (f"{r.get('MAG_Identity', 0):.1f}%"
                 if pd.notna(r.get('MAG_Identity')) and r.get('MAG_Identity', 0) > 0
                 else '-')
        mag_b = (str(r.get('Bin', ''))
                 if pd.notna(r.get('Bin')) else '-')
        s_name = str(r.get('Sample', ''))
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{s_name}">{s_name}</td>'
            f'<td><strong>{_esc(r["Gene"])}</strong></td>'
            f'<td>{badge}</td>'
            f'<td>{kma_d}</td><td>{kma_i}</td>'
            f'<td>{mag_i}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Gene</th>'
        '<th>Detection</th><th>KMA Depth</th><th>KMA Identity</th>'
        '<th>MAG Identity</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Reads+Contigs = maximum confidence. Reads only = possible '
        'unassembled reservoir. Contigs only = insufficient read depth for KMA.</p>')


def html_vfdb_table(vfdb_df):
    if vfdb_df.empty:
        return '<p class="no-data">No virulence factors detected.</p>'
    rows = ''
    for _, r in vfdb_df.iterrows():
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{r["Sample"]}">{r["Sample"]}</td>'
            f'<td>{_esc(r["Bin"])}</td>'
            f'<td><strong>{_esc(r["Gene"])}</strong></td>'
            f'<td class="desc-cell">{_esc(r["Product"])}</td>'
            f'<td>{r["Identity"]:.1f}%</td>'
            f'<td>{r["Coverage"]:.1f}%</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Bin</th><th>Gene</th><th>Product</th>'
        '<th>% Identity</th><th>% Coverage</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        f'<p class="note">VFDB database (Virulence Factor Database). '
        f'{len(vfdb_df)} total hits detected by ABRicate.</p>')


def _card_search_url(query):
    """CARD ARO search URL — works for gene names (otr(B), fos, blaVIM,
    vanR-B …) and returns the list of related ARO entries with their
    full ontology metadata. NCBI Pathogens RefGene was tried first but
    returns unrelated hits when fed gene families like `fos` or `tet`.
    CARD's ARO browser treats short gene names as a family-level search
    and returns the actual variants — closer to what a clinician needs.
    """
    import urllib.parse as up
    return ('https://card.mcmaster.ca/aro/?search_query='
            + up.quote_plus(query))


def _refgene_search_url(query):  # backward-compat shim
    return _card_search_url(query)


def _refgene_url(nf_id):  # backward-compat shim
    return _card_search_url(nf_id.split('.', 1)[0])


def _enrich_genomad_amr_ids(amr_id_str, amr_plasmids_df, contig,
                            amrfinder_catalog=None):
    """geNomad's `amr_genes` column lists NF identifiers (NCBI Reference
    Gene catalogue, e.g. NF033135). Each ID is annotated with the
    canonical AMRFinderPlus catalogue (fam.tsv) when available — gene
    symbol, AMR class, AMR subclass and family/product — and as a
    secondary source we also pull anything detected by the AMRFinderPlus-
    on-plasmids pass on the same contig.

    Returns a small HTML block per NF id with an active RefGene link
    and a tooltip carrying the family/product description.
    """
    if not amr_id_str or amr_id_str.strip() in ('', 'NA'):
        return ''
    ids = [x.strip() for x in re.split(r'[;,]', amr_id_str) if x.strip()]

    # Detected-on-this-contig lookup from AMRFinderPlus-on-plasmids
    contig_lookup = {}
    if amr_plasmids_df is not None and not amr_plasmids_df.empty:
        for _, ar in amr_plasmids_df[amr_plasmids_df['Contig'] == contig].iterrows():
            acc = str(ar.get('Accession', '')).strip().split('.', 1)[0]
            contig_lookup[acc] = {
                'gene_symbol': str(ar.get('Gene', '')).strip(),
                'class'      : str(ar.get('Class', '')).strip(),
                'subclass'   : str(ar.get('Subclass', '')).strip(),
                'product'    : str(ar.get('Name', '')).strip(),
            }

    chips = []
    catalog = amrfinder_catalog or {}
    for nf_id in ids:
        base = nf_id.split('.', 1)[0]
        info = dict(catalog.get(base, {}))   # canonical first
        # Layer in the AMRFinderPlus-on-plasmids hit for the same
        # contig — keeps the catalogue defaults when the on-plasmid pass
        # left a field empty.
        for k, v in contig_lookup.get(base, {}).items():
            if v and not info.get(k):
                info[k] = v

        gene = info.get('gene_symbol', '')
        cls = info.get('class', '')
        sub = info.get('subclass', '')
        family = info.get('family', '') or info.get('product', '')

        cls_chip = ''
        if cls:
            color = _amr_class_color(cls)
            cls_chip = (f'<span class="amr-class-chip" '
                        f'style="background:{color}">{_esc(cls.title())}</span>')
            if sub and sub != cls:
                cls_chip += (f'<span class="amr-class-chip" '
                             f'style="background:#5B6770">'
                             f'{_esc(sub.title())}</span>')

        title_attr = (f'CARD ARO search · {_esc(gene)} · {_esc(nf_id)}'
                      + (f' · {_esc(family)}' if family else ''))
        link_label = (
            f'<strong>{_esc(gene)}</strong>' if gene
            else f'<code>{_esc(nf_id)}</code>')

        chip_html = (
            f'<span class="amr-chip">'
            f'<a href="{_refgene_url(nf_id)}" target="_blank" rel="noopener" '
            f'title="{title_attr}">{link_label}</a>'
            f'{cls_chip}'
            f'</span>')
        chips.append(chip_html)
    return ' '.join(chips)


def html_plasmid_table(plasmids_df, integ_df, amr_plasmids_df=None,
                       amrfinder_catalog=None):
    if plasmids_df.empty:
        return '<p class="no-data">No plasmids detected by geNomad.</p>'

    df = plasmids_df.copy()
    df['Plasmid_score'] = pd.to_numeric(df['Plasmid_score'], errors='coerce').fillna(0)
    df['Length'] = pd.to_numeric(df['Length'], errors='coerce').fillna(0).astype(int)
    df['N_genes'] = pd.to_numeric(df['N_genes'], errors='coerce').fillna(0).astype(int)
    df['AMR_genes'] = df.get('AMR_genes', '').fillna('').astype(str)

    # ── Per-sample summary (always shown) ──
    summary_rows = ''
    for s in sorted(df['Sample'].unique()):
        s_df = df[df['Sample'] == s]
        n_plas = len(s_df)
        n_conf = int((s_df['Plasmid_score'] >= 0.9).sum())
        with_amr = s_df[s_df['AMR_genes'].str.strip().ne('')
                        & s_df['AMR_genes'].str.strip().ne('NA')]
        n_amr_plas = len(with_amr)
        summary_rows += (
            f'<tr><td class="sample" data-sample="{s}">{s}</td>'
            f'<td>{n_plas:,}</td><td>{n_conf:,}</td><td>{n_amr_plas:,}</td></tr>')
    summary_table = (
        '<h4>Per-sample summary</h4>'
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Plasmid contigs (total)</th>'
        '<th>Confirmed (score &ge;0.9)</th>'
        '<th>Carrying AMR (geNomad)</th>'
        f'</tr></thead><tbody>{summary_rows}</tbody></table>')

    # ── Main detail: ONLY plasmid contigs that carry AMR ──
    amr_only = df[df['AMR_genes'].str.strip().ne('')
                  & df['AMR_genes'].str.strip().ne('NA')].copy()
    amr_only = amr_only.sort_values(['Sample', 'Plasmid_score'],
                                    ascending=[True, False])
    if amr_only.empty:
        amr_only_table = ('<p class="note">No plasmid contig carries AMR '
                          'annotated by geNomad in this run. The independent '
                          'AMRFinderPlus pass below may still find hits.</p>')
    else:
        # Build one row per (contig, NF id) with explicit columns. The
        # canonical info comes from the AMRFinderPlus reference catalogue
        # (fam.tsv) — same data NCBI publishes for that DB version. The
        # AMRFinderPlus-on-plasmids pass is layered in only to fill in
        # any field the catalogue left empty for that contig.
        contig_lookup = {}
        if amr_plasmids_df is not None and not amr_plasmids_df.empty:
            for _, ar in amr_plasmids_df.iterrows():
                acc = str(ar.get('Accession', '')).strip().split('.', 1)[0]
                contig_lookup[(str(ar.get('Contig', '')), acc)] = {
                    'gene_symbol': str(ar.get('Gene', '')).strip(),
                    'class'      : str(ar.get('Class', '')).strip(),
                    'subclass'   : str(ar.get('Subclass', '')).strip(),
                    'product'    : str(ar.get('Name', '')).strip(),
                }
        catalog = amrfinder_catalog or {}

        rows = ''
        n_lines = 0
        for _, r in amr_only.iterrows():
            score = float(r['Plasmid_score'])
            score_cls = ('good' if score >= 0.9
                         else ('warn' if score >= 0.7 else 'bad'))
            ids = [x.strip() for x in re.split(r'[;,]', r['AMR_genes']) if x.strip()]
            for nf_id in ids:
                base = nf_id.split('.', 1)[0]
                info = dict(catalog.get(base, {}))
                for k, v in contig_lookup.get((r['Contig'], base), {}).items():
                    if v and not info.get(k):
                        info[k] = v
                gene = info.get('gene_symbol', '') or '<em>—</em>'
                cls = info.get('class', '')
                sub = info.get('subclass', '')
                family = info.get('family', '') or info.get('product', '')

                cls_chip = ''
                if cls:
                    color = _amr_class_color(cls)
                    cls_chip = (f'<span class="amr-class-chip" '
                                f'style="background:{color}">'
                                f'{_esc(cls.title())}</span>')
                sub_cell = (
                    f'<small>{_esc(sub.replace("/", " / "))}</small>'
                    if sub and sub != cls else '')

                gene_clean = info.get('gene_symbol', '')
                if gene_clean:
                    gene_cell = (
                        f'<a href="{_refgene_search_url(gene_clean)}" '
                        f'target="_blank" rel="noopener" '
                        f'title="CARD ARO search for {_esc(gene_clean)} variants">'
                        f'<strong>{_esc(gene_clean)}</strong></a>')
                else:
                    gene_cell = '<em>—</em>'

                rows += (
                    f'<tr><td class="sample" data-sample="{r["Sample"]}">'
                    f'{_esc(r["Sample"])}</td>'
                    f'<td><code>{_esc(r["Contig"])}</code></td>'
                    f'<td>{int(r["Length"]):,}</td>'
                    f'<td class="{score_cls}">{score:.3f}</td>'
                    f'<td>{gene_cell}</td>'
                    f'<td>{cls_chip}</td>'
                    f'<td>{sub_cell}</td>'
                    f'<td class="desc-cell">{_esc(family)}</td>'
                    f'<td><code title="Internal AMRFinderPlus identifier — '
                    f'NCBI does not expose a stable URL for individual '
                    f'NF ids; use the Gene column to look up the family.">'
                    f'{_esc(nf_id)}</code></td></tr>')
                n_lines += 1

        amr_only_table = (
            f'<h4>Plasmid contigs carrying AMR ({len(amr_only)} contigs · '
            f'{n_lines} gene calls)</h4>'
            '<table class="data-table sortable filterable"><thead><tr>'
            '<th>Sample</th><th>Contig</th><th>Length (bp)</th>'
            '<th>Plasmid score</th>'
            '<th>Gene</th><th>Class</th><th>Subclass</th>'
            '<th>Family / product</th><th>NF id</th>'
            '</tr></thead>'
            f'<tbody>{rows}</tbody></table>'
            '<p class="note"><strong>Gene · Class · Subclass · Family/product</strong> '
            'come from the AMRFinderPlus reference catalogue '
            '(<code>fam.tsv</code>) shipped with the NCBI database for '
            'this run — the canonical source. '
            '<strong>The Gene column links to a CARD ARO search</strong> by '
            'gene name, which returns the family with all its variants and '
            'their ontology — useful when you need to confirm whether a '
            'short family name (<code>fos</code>, <code>tet</code>, '
            '<code>bla</code>) matches a specific clinically relevant '
            'allele. NCBI Pathogens RefGene was tried first but returned '
            'unrelated records when fed family-level names. '
            '<strong>The NF id is shown as a code without a link</strong>: '
            'no public catalogue exposes a stable URL for individual '
            'NF identifiers, so the catalogue columns above are the '
            'canonical reference and the link is a starting point. '
            'One row per NF id per contig.</p>')

    # ── All plasmid contigs (potentially thousands) — drill-down only ──
    all_rows = ''
    head_df = df.sort_values(['Sample', 'Plasmid_score'],
                             ascending=[True, False]).head(2000)
    for _, r in head_df.iterrows():
        score = float(r['Plasmid_score'])
        score_cls = ('good' if score >= 0.9
                     else ('warn' if score >= 0.7 else 'bad'))
        amr_str = str(r.get('AMR_genes', '')) or '-'
        all_rows += (
            f'<tr><td class="sample" data-sample="{r["Sample"]}">'
            f'{_esc(r["Sample"])}</td>'
            f'<td><code>{_esc(r["Contig"])}</code></td>'
            f'<td>{int(r["Length"]):,}</td>'
            f'<td class="{score_cls}">{score:.3f}</td>'
            f'<td>{int(r["N_genes"])}</td>'
            f'<td>{_esc(amr_str)}</td></tr>')
    truncation_note = ''
    if len(df) > 2000:
        truncation_note = (
            f'<p class="note">First 2,000 plasmid contigs of {len(df):,} '
            'total are shown. The full list lives in '
            '<code>26_genomad/&lt;sample&gt;/&lt;sample&gt;_plasmid_summary.tsv</code>.</p>')
    all_table = (
        '<details class="prov-panel" style="margin-top:0.6rem">'
        '<summary><span class="prov-icon">+</span>'
        '<span>Drill-down: all plasmid contigs predicted by geNomad '
        f'(showing up to 2,000 of {len(df):,})</span></summary>'
        + truncation_note +
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Contig</th><th>Length (bp)</th>'
        '<th>Plasmid score</th><th>N genes</th>'
        '<th>AMR genes (raw NF ids)</th>'
        f'</tr></thead><tbody>{all_rows}</tbody></table></details>')

    return (
        summary_table + amr_only_table + all_table +
        '<p class="note">Plasmid score (geNomad): &ge;0.9 confirmed, '
        '0.7–0.9 probable, &lt;0.7 uncertain. The full TSV is the canonical '
        'artefact for downstream analysis — the table above filters to the '
        'rows that carry actionable AMR signal.</p>')


def html_integron_table(integron_summary_df, amr_mags_df=None, tax_df=None):
    """Per-MAG integron summary with cross-link to AMR carriage and
    organism. Honest when no integrons were found in any bin."""
    if integron_summary_df is None or integron_summary_df.empty:
        return ('<p class="no-data"><strong>IntegronFinder did not detect '
                'any integron in any of the recovered MAGs for this run.</strong> '
                'This is consistent with bins of organisms (e.g. '
                '<em>Acinetobacter johnsonii</em>) that lack the canonical '
                'class-1/2/3 integron platform, or with fragmented bins '
                'where the integrase + cassette layout is broken across '
                'contigs. The integrase-aware <code>intI</code> annotation '
                'is the highest-confidence signal — its absence here is '
                'biological, not a bug.</p>')

    # Cross-links: AMR genes per (Sample, Bin) and organism per (Sample, Bin)
    amr_by_bin = {}
    if amr_mags_df is not None and not amr_mags_df.empty:
        for _, r in amr_mags_df.iterrows():
            key = (str(r.get('Sample', '')), str(r.get('Bin', '')))
            amr_by_bin.setdefault(key, []).append(
                f'{r.get("Gene", "")} ({r.get("Class", "")})')
    organism_by_bin = {}
    if tax_df is not None and not tax_df.empty:
        for _, r in tax_df.iterrows():
            key = (str(r.get('Sample', '')), str(r.get('Bin', '')))
            organism_by_bin[key] = str(r.get('Organism', ''))

    rows = ''
    for _, r in integron_summary_df.iterrows():
        has_intI_badge = (
            f'<span class="badge" style="background:{RED}">Yes</span>'
            if r['has_intI'] == 'Yes'
            else f'<span class="badge" style="background:{GREEN}">No</span>')
        bin_key = (str(r['Sample']), str(r['Bin']))
        organism = organism_by_bin.get(bin_key, '')
        amr_in_bin = sorted(set(amr_by_bin.get(bin_key, [])))
        amr_html = (', '.join(amr_in_bin) if amr_in_bin
                    else '<span class="note">none in this bin</span>')
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{r["Sample"]}">{r["Sample"]}</td>'
            f'<td>{_esc(r["Bin"])}</td>'
            f'<td><em>{_esc(organism)}</em></td>'
            f'<td>{r["N_complete"]}</td>'
            f'<td>{r["N_In0"]}</td>'
            f'<td>{r["N_CALIN"]}</td>'
            f'<td>{has_intI_badge}</td>'
            f'<td>{amr_html}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>MAG</th><th>Organism</th>'
        '<th>Complete</th><th>In0</th><th>CALIN</th>'
        '<th>Integrase (intI)</th><th>AMR genes in this bin</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Complete = integrase + cassettes (active platform); '
        'In0 = integrase only (capable but empty); CALIN = cassettes without '
        'integrase (orphan). The AMR-genes column is a cross-link with '
        'AMRFinderPlus on the same MAG — co-occurrence of <code>intI</code> '
        'and AMR genes in a bin is the signal of an active mobilisation '
        'platform, even when IntegronFinder cannot place specific cassettes.</p>')


_TAX_TOOL_KEY = {
    'Skani'   : 'skani',
    'Sourmash': 'sourmash',
    'GTDB-Tk' : 'gtdbtk',
}

# Pipeline-order display so the Software & Versions table reads as a flow
SOFTWARE_TABLE_ORDER = [
    'nanoplot', 'fastqc', 'fastqscreen', 'porechop', 'chopper',
    'minimap2_host', 'samtools', 'kraken2', 'bracken', 'kaiju',
    'sylph', 'kma', 'minimap2_phenotypic', 'multiqc',
    'metaflye', 'medaka', 'quast',
    'metabat2', 'maxbin2', 'semibin2', 'dastool',
    'checkm2', None,  # placeholder for the active genome-taxonomy tool
    'amrfinderplus', 'genomad', 'abricate', 'mobsuite',
    'integronfinder', 'bakta',
]


def html_software_table(tax_info):
    """Software & versions table — pipeline-ordered, picks the active
    genome-taxonomy classifier from `tax_info`."""
    active_tax_key = _TAX_TOOL_KEY.get(tax_info['label'], 'skani')
    rows = ''
    for key in SOFTWARE_TABLE_ORDER:
        if key is None:
            # Slot for the genome-taxonomy classifier, with DB version inline
            label, version, role, _layer, url, _, _ = TOOL_REGISTRY[active_tax_key]
            version_full = f'{version} ({tax_info["db_label"]})'
            rows += (
                f'<tr><td><strong>{_esc(label)}</strong></td>'
                f'<td>{_esc(version_full)}</td>'
                f'<td>{_esc(role)}</td>'
                f'<td><a href="{url}" target="_blank">{_esc(url)}</a></td></tr>\n')
            continue
        if key not in TOOL_REGISTRY:
            continue
        label, version, role, _layer, url, _, _ = TOOL_REGISTRY[key]
        rows += (
            f'<tr><td><strong>{_esc(label)}</strong></td>'
            f'<td>{_esc(version)}</td>'
            f'<td>{_esc(role)}</td>'
            f'<td><a href="{url}" target="_blank">{_esc(url)}</a></td></tr>\n')
    return (
        '<table class="data-table sortable"><thead><tr>'
        '<th>Software</th><th>Version</th><th>Function</th><th>URL</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>')


def _format_param_value(v):
    if v is None:
        return '<em>not set</em>'
    s = str(v)
    if len(s) > 90:
        return f'<code title="{_esc(s)}">{_esc(s[:80])}…</code>'
    return f'<code>{_esc(s)}</code>'


def _runtime_summary_for_tool(trace_df, tool_key):
    """Compact runtime string for a tool, derived from trace.txt rows.
    For multi-task tools (Bakta per-bin, NanoPlot raw+filtered, etc.) we
    show n executions plus min–max range and mean realtime instead of
    a single sum, so per-sample variability is visible."""
    sub = trace_for_tool(trace_df, tool_key)
    if sub.empty:
        return ''
    n = len(sub)
    n_cached = int((sub['status'] == 'CACHED').sum())
    total_realtime = float(sub['Realtime_s'].sum())
    peak = float(sub['Peak_rss_b'].max())
    cached_note = (f'{n_cached} cached / {n - n_cached} run' if n_cached else
                   f'{n} run')

    if n == 1:
        runtime_str = _fmt_seconds(total_realtime)
    else:
        rt_min = float(sub['Realtime_s'].min())
        rt_max = float(sub['Realtime_s'].max())
        rt_mean = total_realtime / n
        runtime_str = (
            f'{_fmt_seconds(total_realtime)} total · '
            f'{n}× · {_fmt_seconds(rt_min)}–{_fmt_seconds(rt_max)} '
            f'(μ {_fmt_seconds(rt_mean)})')

    return (f'<a href="#sec-execution" title="Open Pipeline execution section">'
            f'⏱ {runtime_str}</a> &middot; '
            f'<span title="Peak RSS across all tasks">'
            f'{_fmt_bytes(peak)} peak</span> &middot; '
            f'<span class="trace-status">{cached_note}</span>')


def html_provenance_panel(section_key, run_params, trace_df=None):
    """Render a collapsible "Source / tools / parameters / runtime" panel
    for a given section. `run_params` is the dict from `load_run_params()`;
    `trace_df` is `load_trace()` output. Returns '' if the section has no
    registered tools."""
    tools = SECTION_TOOLS.get(section_key, [])
    if not tools:
        return ''
    params_dict = (run_params or {}).get('params', {}) or {}
    if trace_df is None:
        trace_df = pd.DataFrame()
    has_runtime = not trace_df.empty

    rows = ''
    for tool_key in tools:
        if tool_key not in TOOL_REGISTRY:
            continue
        label, version, role, layer, url, param_keys, _proc = TOOL_REGISTRY[tool_key]
        param_lines = []
        for pk in param_keys:
            if pk in params_dict:
                param_lines.append(
                    f'<li><code>{_esc(pk)}</code>: {_format_param_value(params_dict[pk])}</li>')
            else:
                param_lines.append(
                    f'<li><code>{_esc(pk)}</code>: <em>default</em></li>')
        params_html = (f'<ul class="prov-params">{"".join(param_lines)}</ul>'
                       if param_lines else '<em>no run-tunable parameters</em>')

        runtime_html = (_runtime_summary_for_tool(trace_df, tool_key)
                        if has_runtime else '')
        runtime_cell = (f'<td class="prov-runtime">{runtime_html}</td>'
                        if has_runtime else '')

        rows += (
            f'<tr><td><a href="{url}" target="_blank" rel="noopener">'
            f'<strong>{_esc(label)}</strong></a></td>'
            f'<td>{_esc(version)}</td>'
            f'<td>{_esc(LAYER_LABELS.get(layer, layer))}</td>'
            f'<td>{_esc(role)}</td>'
            f'<td>{params_html}</td>'
            f'{runtime_cell}</tr>\n')

    runtime_th = '<th>Runtime (this run)</th>' if has_runtime else ''
    return (
        '<details class="prov-panel"><summary>'
        '<span class="prov-icon">i</span> '
        '<span>Source &amp; parameters for this section</span>'
        '</summary>'
        '<table class="prov-table"><thead><tr>'
        '<th>Tool</th><th>Version</th><th>Layer</th><th>Role</th>'
        f'<th>Parameters used</th>{runtime_th}</tr></thead>'
        f'<tbody>{rows}</tbody></table></details>')


def html_trace_summary_cards(summary, run_params):
    """Top cards for the Pipeline-execution section."""
    if not summary:
        return '<p class="no-data">No trace.txt found for this run. ' \
               'Re-run the pipeline so Nextflow writes the trace.</p>'
    duration = (run_params or {}).get('run', {}).get('duration', '—')
    cards = [
        ('Wall-clock duration', _esc(duration), 'Total run time (start → end)'),
        ('Total processes',     str(summary['n_total']),
         f'{summary["n_completed"]} executed · {summary["n_cached"]} cached'),
        ('Reused (cached)',     f'{summary["pct_cached"]:.0f}%',
         'Higher = better incremental rebuild'),
        ('CPU-hours',           f'{summary["cpu_hours"]:.1f}',
         'Sum of realtime across all tasks'),
        ('Failed tasks',        str(summary['n_failed']),
         'Should be zero on a successful run'),
    ]
    items = ''.join(
        f'<div class="trace-card">'
        f'<div class="trace-card-value">{value}</div>'
        f'<div class="trace-card-label">{_esc(label)}</div>'
        f'<div class="trace-card-hint">{_esc(hint)}</div>'
        f'</div>'
        for label, value, hint in cards)
    return f'<div class="trace-cards">{items}</div>'


def html_trace_top_table(trace_df, by, n, title):
    df = trace_top_n(trace_df, by, n=n,
                     status_filter=['COMPLETED'])
    if df.empty:
        return f'<h4>{_esc(title)}</h4><p class="no-data">No completed tasks.</p>'
    rows = ''
    for _, r in df.iterrows():
        rows += (
            f'<tr><td><code>{_esc(r["Process_short"])}</code></td>'
            f'<td>{_esc(str(r.get("tag", "")))}</td>'
            f'<td>{_esc(str(r["status"]))}</td>'
            f'<td>{_fmt_seconds(r["Realtime_s"])}</td>'
            f'<td>{_fmt_bytes(r["Peak_rss_b"])}</td>'
            f'<td>{_esc(str(r.get("%cpu", "")))}</td></tr>')
    return (
        f'<h4>{_esc(title)}</h4>'
        '<table class="prov-table"><thead><tr>'
        '<th>Process</th><th>Tag</th><th>Status</th>'
        '<th>Realtime</th><th>Peak RSS</th><th>%CPU</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>')


_ANOMALY_BADGE = {
    'RAM'      : f'<span class="badge" style="background:{ORANGE}">RAM</span>',
    'SLOW'     : f'<span class="badge" style="background:{RED}">SLOW</span>',
    'RAM-P95'  : f'<span class="badge" style="background:#9B59B6">RAM-P95</span>',
    'SLOW-P95' : f'<span class="badge" style="background:#3498DB">SLOW-P95</span>',
}


def html_trace_anomalies(trace_df):
    df = trace_anomalies(trace_df)
    if df.empty:
        return ('<p class="no-data">No tasks were flagged as anomalous '
                f'(thresholds: &gt;{TRACE_ANOMALY_RAM_GB} GB peak RSS, '
                f'&gt;{TRACE_ANOMALY_DURATION_S // 60} min realtime, '
                'or &gt;P95 within this run).</p>')
    df = df.sort_values(['Realtime_s', 'Peak_rss_b'], ascending=False)
    rows = ''
    for _, r in df.iterrows():
        flags = ' '.join(_ANOMALY_BADGE.get(f, f) for f in r['Anomaly_flags'])
        rows += (
            f'<tr><td><code>{_esc(r["Process_short"])}</code></td>'
            f'<td>{_esc(str(r.get("tag", "")))}</td>'
            f'<td>{flags}</td>'
            f'<td>{_fmt_seconds(r["Realtime_s"])}</td>'
            f'<td>{_fmt_bytes(r["Peak_rss_b"])}</td></tr>')
    return (
        '<table class="prov-table"><thead><tr>'
        '<th>Process</th><th>Tag</th><th>Flags</th>'
        '<th>Realtime</th><th>Peak RSS</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>'
        '<p class="note"><strong>RAM</strong> = peak RSS above '
        f'{TRACE_ANOMALY_RAM_GB} GB · '
        f'<strong>SLOW</strong> = realtime above {TRACE_ANOMALY_DURATION_S // 60} min · '
        '<strong>RAM-P95 / SLOW-P95</strong> = within the top 5% of '
        'this run (relative outlier, regardless of absolute level).</p>')


def chart_trace_distribution(trace_df, metric, title, axis_title,
                             log_scale=False, value_formatter=None):
    """Distribution chart per process. Violin per process when there are
    enough tasks; falls back to a strip plot when too few. Log scale for
    realtime/RAM as requested.
    `metric` is the dataframe column with numeric values."""
    if trace_df.empty:
        return None
    df = trace_df.copy()
    df = df[df[metric] > 0]
    if df.empty:
        return None

    # Pick the top-12 processes by count to keep the chart legible.
    counts = df['Process_short'].value_counts()
    keep = counts.head(12).index.tolist()
    df = df[df['Process_short'].isin(keep)]

    # Order processes by their median value for the metric, descending.
    medians = df.groupby('Process_short')[metric].median().sort_values(
        ascending=True)
    proc_order = medians.index.tolist()

    fig = go.Figure()
    layer_color = {'TAX': BLUE, 'MAG': ORANGE, 'SETUP': GREEN}
    for proc in proc_order:
        sub = df[df['Process_short'] == proc]
        layer = sub['Process_layer'].iloc[0] if not sub.empty else 'OTHER'
        color = layer_color.get(layer, GREY)
        # Use a violin only if at least 3 points, otherwise strip
        n = len(sub)
        hover = sub.apply(
            lambda r: (f'{r["Process_short"]} ({r["Process_layer"]})<br>'
                       f'tag={r.get("tag", "-")}<br>'
                       f'status={r["status"]}<br>'
                       f'realtime={_fmt_seconds(r["Realtime_s"])}<br>'
                       f'peak RSS={_fmt_bytes(r["Peak_rss_b"])}<br>'
                       f'%cpu={r.get("%cpu", "-")}'), axis=1)
        if n >= 3:
            fig.add_trace(go.Violin(
                y=[proc] * n, x=sub[metric], name=proc, orientation='h',
                marker=dict(color=color), points='all',
                hovertext=hover, hoverinfo='text',
                line=dict(color=color), box_visible=True,
                meanline_visible=True))
        else:
            fig.add_trace(go.Scatter(
                y=[proc] * n, x=sub[metric], mode='markers',
                marker=dict(color=color, size=10, line=dict(color='white', width=1)),
                name=proc, hovertext=hover, hoverinfo='text'))

    fig.update_layout(
        title=dict(text=title, font_color=BLUE, font_size=14),
        height=max(280, 26 * len(proc_order) + 100),
        xaxis_title=axis_title,
        margin=dict(l=200, r=30, t=60, b=60),
        showlegend=False,
        font=dict(family='Arial, sans-serif', size=10),
    )
    if log_scale:
        fig.update_xaxes(type='log')
    return fig


def html_trace_distributions(trace_df):
    """Three resource-distribution charts arranged vertically. Empty
    when trace.txt is unavailable."""
    if trace_df.empty:
        return ''
    rt_chart = fig_to_html(chart_trace_distribution(
        trace_df, 'Realtime_s',
        'Realtime distribution (log scale)',
        'Realtime (s, log)', log_scale=True))
    ram_chart = fig_to_html(chart_trace_distribution(
        trace_df, 'Peak_rss_b',
        'Peak RSS distribution (log scale)',
        'Peak RSS (bytes, log)', log_scale=True))
    # %cpu stays linear; values are bounded ~0–CPUs*100
    cpu_df = trace_df.copy()
    cpu_df['_cpu_num'] = pd.to_numeric(
        cpu_df.get('%cpu', '').astype(str).str.replace('%', '', regex=False),
        errors='coerce').fillna(0)
    cpu_chart = ''
    if not cpu_df.empty and cpu_df['_cpu_num'].sum() > 0:
        # Reuse the chart helper by swapping the metric name in.
        cpu_df['Realtime_s'] = cpu_df['_cpu_num']
        cpu_chart = fig_to_html(chart_trace_distribution(
            cpu_df, 'Realtime_s',
            '%CPU distribution (linear)',
            '%CPU (single-CPU = 100)', log_scale=False))
    return (rt_chart + ram_chart + cpu_chart) if (rt_chart or ram_chart) else ''


def html_trace_by_layer(trace_df):
    """Group trace rows by layer (SETUP / TAX / MAG) and render as
    sortable tables. Each row has process, tag, status, realtime, peak."""
    if trace_df.empty:
        return ''
    blocks = []
    for layer in TRACE_LAYERS:
        sub = trace_df[trace_df['Process_layer'] == layer].sort_values(
            'Realtime_s', ascending=False)
        if sub.empty:
            continue
        sub_n_run = int((sub['status'] == 'COMPLETED').sum())
        sub_n_cached = int((sub['status'] == 'CACHED').sum())
        rows = ''
        for _, r in sub.iterrows():
            status_cls = ('good' if r['status'] == 'COMPLETED'
                          else 'note')
            rows += (
                f'<tr><td><code>{_esc(r["Process_short"])}</code></td>'
                f'<td>{_esc(str(r.get("tag", "")))}</td>'
                f'<td class="{status_cls}">{_esc(str(r["status"]))}</td>'
                f'<td>{_fmt_seconds(r["Realtime_s"])}</td>'
                f'<td>{_fmt_bytes(r["Peak_rss_b"])}</td>'
                f'<td>{_esc(str(r.get("%cpu", "")))}</td></tr>')
        blocks.append(
            f'<h4>{_esc(layer)} layer · {len(sub)} tasks '
            f'({sub_n_cached} cached, {sub_n_run} run)</h4>'
            '<table class="data-table sortable filterable">'
            '<thead><tr>'
            '<th>Process</th><th>Tag</th><th>Status</th>'
            '<th>Realtime</th><th>Peak RSS</th><th>%CPU</th>'
            f'</tr></thead><tbody>{rows}</tbody></table>')
    return ''.join(blocks)


def html_run_parameters_panel(run_params):
    """Top-level run-parameters table: the full params.json contents,
    grouped by category for readability. Falls back to a notice if
    params.json was not produced (legacy run)."""
    if not run_params or 'params' not in run_params:
        return ('<p class="no-data">No params.json was found for this '
                'run. Re-run the pipeline with the current main.nf to '
                'enable per-run parameter trace.</p>')

    run_meta = run_params.get('run', {}) or {}
    params = run_params.get('params', {}) or {}

    meta_rows = ''
    for k in ['name', 'started_at', 'completed_at', 'duration', 'success',
              'profile', 'nextflow', 'revision', 'session_id', 'command_line']:
        if k in run_meta and run_meta[k] is not None:
            meta_rows += (f'<tr><td><code>{_esc(k)}</code></td>'
                          f'<td>{_format_param_value(run_meta[k])}</td></tr>')
    meta_table = (
        '<h4>Run metadata</h4>'
        '<table class="prov-table"><thead><tr><th>Field</th><th>Value</th>'
        f'</tr></thead><tbody>{meta_rows}</tbody></table>') if meta_rows else ''

    # Categorise pipeline parameters
    categories = [
        ('Inputs',     ['input', 'outdir', 'run_name', 'run_tax', 'run_assembly',
                        'taxonomy_tool', 'filtered_dir', 'assemblies_dir']),
        ('Filtering',  ['min_quality', 'min_length', 'host_genome',
                        'phenotypic_targets', 'min_coverage']),
        ('Databases',  ['db_root', 'kraken2_db', 'kaiju_nodes', 'kaiju_names',
                        'resfinder_fsa', 'fastqscreen_conf', 'sylph_db',
                        'checkm2_db', 'gtdbtk_data', 'bakta_db', 'bakta_db_type',
                        'genomad_db', 'amrfinder_db', 'skani_db']),
        ('Other',      []),
    ]
    seen = set()
    for _, keys in categories:
        for k in keys:
            seen.add(k)
    others = sorted(k for k in params if k not in seen)
    categories[-1] = ('Other', others)

    cat_html = ''
    for cat_name, keys in categories:
        rows = ''
        for k in keys:
            if k in params:
                rows += (f'<tr><td><code>{_esc(k)}</code></td>'
                         f'<td>{_format_param_value(params[k])}</td></tr>')
        if not rows:
            continue
        cat_html += (
            f'<h4>{_esc(cat_name)}</h4>'
            '<table class="prov-table"><thead><tr><th>Parameter</th>'
            f'<th>Value</th></tr></thead><tbody>{rows}</tbody></table>')

    return meta_table + cat_html


# ═══════════════════════════════════════════════════════════════
# CSS
# ═══════════════════════════════════════════════════════════════

HTML_CSS = """<style>
:root {
    --blue: #2C5F8A; --orange: #E05C2A; --green: #28A745;
    --red: #E74C3C; --grey: #6C757D; --bg: #F8F9FA;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
    font-family: 'Segoe UI', Tahoma, Arial, sans-serif;
    background: var(--bg); color: #333; line-height: 1.6;
}
.header {
    background: linear-gradient(135deg, var(--blue) 0%, #1a3f5c 100%);
    color: white; padding: 1.5rem 2rem;
    border-bottom: 4px solid var(--orange);
    display: flex; align-items: center; gap: 1.5rem;
}
.header-logo { max-height: 45px; width: auto; }
.header-text h1 { font-size: 1.5rem; margin-bottom: 0.2rem; }
.header-text .subtitle { opacity: 0.8; font-size: 0.85rem; }
.header-text .org-info { opacity: 0.6; font-size: 0.75rem; }
.container { max-width: 1400px; margin: 0 auto; padding: 1.5rem; }
.metrics-row {
    display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 1.5rem;
}
.metric-card {
    flex: 1; min-width: 130px; background: white; border-radius: 8px;
    padding: 1rem; box-shadow: 0 2px 6px rgba(0,0,0,0.08);
    text-align: center; border-top: 3px solid var(--blue);
}
.metric-card .value {
    font-size: 1.8rem; font-weight: 700; color: var(--blue);
}
.metric-card .label {
    font-size: 0.78rem; color: var(--grey);
}
.section {
    background: white; border-radius: 8px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    margin-bottom: 1.5rem; overflow: hidden;
}
.section-header {
    background: var(--blue); color: white;
    padding: 0.8rem 1.5rem; font-size: 1.05rem; font-weight: 600;
}
.section-body { padding: 1.5rem; overflow-x: auto; }
.section-body > p {
    margin-bottom: 0.8rem; color: #555; font-size: 0.9rem;
}
.data-table {
    width: 100%; border-collapse: collapse;
    font-size: 0.82rem; margin: 1rem 0;
}
.data-table th {
    background: var(--blue); color: white;
    padding: 0.5rem 0.6rem; text-align: center;
    font-weight: 600; font-size: 0.78rem;
    cursor: pointer; position: relative; user-select: none;
    white-space: nowrap;
}
.data-table th:hover { background: #1a3f5c; }
.data-table td {
    padding: 0.4rem 0.6rem; text-align: center;
    border-bottom: 1px solid #eee;
}
.data-table tr:nth-child(even) { background: #F0F4F7; }
.data-table tr:hover { background: #E8EFF5; }
.data-table .sample {
    text-align: left; font-weight: 600;
    font-family: monospace; font-size: 0.78rem;
}
.data-table .good { color: var(--green); font-weight: 600; }
.data-table .warn { color: var(--orange); font-weight: 600; }
.data-table .bad { color: var(--red); font-weight: 600; }
.data-table .highlight {
    background: #C8E6C9 !important; font-weight: 700;
}
.data-table .total-row {
    background: #E3EBF3 !important; font-size: 0.85rem;
}
.data-table .desc-cell {
    text-align: left; font-size: 0.76rem;
    max-width: 250px; word-wrap: break-word;
}
.data-table .links-cell {
    font-size: 0.74rem; white-space: nowrap;
}
.data-table .links-cell a {
    color: var(--blue); text-decoration: none;
}
.data-table .links-cell a:hover { text-decoration: underline; }
.badge {
    display: inline-block; padding: 0.15rem 0.5rem;
    border-radius: 10px; font-size: 0.7rem;
    font-weight: 700; color: white;
}
.note { font-size: 0.78rem; color: var(--grey); margin-top: 0.8rem; }
.no-data {
    color: #999; font-style: italic; padding: 1rem 0;
}
.alert-box {
    background: #FFEBEE; border-left: 4px solid var(--red);
    padding: 1rem; border-radius: 4px; margin: 1rem 0;
}
.alert-box p { margin-bottom: 0.5rem; }
.conclusion {
    background: #E8F5E9; border-left: 4px solid var(--green);
    padding: 1rem; border-radius: 4px;
}
.footer {
    text-align: center; padding: 1.5rem;
    color: var(--grey); font-size: 0.75rem;
}
.sample-filter {
    background: #f0f4f7; padding: 0.8rem 1.2rem;
    border-radius: 6px; margin-bottom: 1rem;
}
.sample-filter label {
    margin-right: 1rem; font-size: 0.85rem; cursor: pointer;
}
.sample-filter input[type=checkbox] { margin-right: 0.3rem; }
.search-box { margin-bottom: 1rem; }
.search-box input {
    width: 100%; padding: 0.5rem 0.8rem;
    border: 1px solid #ccc; border-radius: 4px;
    font-size: 0.85rem;
}
.search-box input:focus {
    outline: none; border-color: var(--blue);
    box-shadow: 0 0 0 2px rgba(44,95,138,0.2);
}
.sort-asc::after { content: ' \\25B2' !important; opacity: 1 !important; }
.sort-desc::after { content: ' \\25BC' !important; opacity: 1 !important; }
.data-table.sortable th::after {
    content: ' \\2195'; font-size: 0.65rem; opacity: 0.5;
}
.dag-container { padding: 1rem; overflow-x: auto; }
.dag-flow {
    display: flex; flex-direction: column;
    gap: 0.3rem; align-items: center; font-size: 0.75rem;
}
.dag-row {
    display: flex; gap: 0.5rem; align-items: center;
    flex-wrap: wrap; justify-content: center;
}
.dag-box {
    padding: 0.4rem 0.7rem; border-radius: 4px;
    color: white; font-weight: 600; white-space: nowrap;
    text-align: center; font-size: 0.72rem;
}
.dag-box.qc { background: #5B9ABF; }
.dag-box.trim { background: #2C5F8A; }
.dag-box.tax { background: #28A745; }
.dag-box.asm { background: #E05C2A; }
.dag-box.bin { background: #9B59B6; }
.dag-box.ann { background: #E74C3C; }
.dag-box.rep { background: #1ABC9C; }
.dag-box.rep-final {
    background: linear-gradient(135deg, #1ABC9C, #2C5F8A);
    font-size: 0.8rem; padding: 0.5rem 1rem;
}
.dag-arrow { color: #aaa; font-size: 1.2rem; }
.dag-arrow-h { color: #aaa; font-size: 0.9rem; }
.dag-label { font-size: 0.65rem; color: #999; margin-top: -0.2rem; }
.dag-step-num {
    font-size: 0.6rem; color: #bbb; font-weight: 600;
    margin-top: 0.3rem;
}
.toc { margin: 1rem 0; }
.toc a {
    color: var(--blue); text-decoration: none;
    display: block; padding: 0.25rem 0 0.25rem 0.5rem;
    font-size: 0.88rem; border-left: 2px solid transparent;
}
.toc a:hover { text-decoration: underline; border-left-color: var(--orange); }
.toc-group-header {
    color: var(--grey); font-size: 0.72rem; font-weight: 700;
    text-transform: uppercase; letter-spacing: 0.05em;
    margin: 0.9rem 0 0.2rem 0;
    padding-bottom: 0.15rem; border-bottom: 1px solid #e0e0e0;
}
.layer-header {
    color: var(--blue); font-size: 1.1rem; font-weight: 600;
    margin: 2.2rem 0 0.8rem 0;
    padding-bottom: 0.4rem;
    border-bottom: 2px solid var(--blue);
    letter-spacing: 0.02em;
}
.prov-panel {
    margin-top: 1rem; background: #f6f8fa;
    border: 1px solid #d8dde2; border-radius: 4px;
    padding: 0.6rem 0.9rem; font-size: 0.85rem;
}
.prov-panel summary {
    cursor: pointer; color: var(--blue); font-weight: 600;
    list-style: none;
}
.prov-panel summary::-webkit-details-marker { display: none; }
.prov-panel[open] summary { margin-bottom: 0.6rem; }
.prov-icon {
    display: inline-block; width: 18px; height: 18px;
    line-height: 18px; text-align: center;
    background: var(--blue); color: white; border-radius: 50%;
    font-style: italic; font-weight: 700; font-size: 0.75rem;
    margin-right: 0.4rem;
}
.prov-table {
    width: 100%; font-size: 0.82rem; border-collapse: collapse;
    margin: 0.4rem 0;
}
.prov-table th {
    background: var(--blue); color: white;
    padding: 0.35rem 0.5rem; text-align: left; font-weight: 600;
}
.prov-table td {
    padding: 0.35rem 0.5rem; border-bottom: 1px solid #e6e6e6;
    vertical-align: top;
}
.prov-params { margin: 0; padding-left: 1rem; font-size: 0.78rem; }
.prov-params li { margin: 0.1rem 0; }
.risk-rule {
    color: var(--grey); font-size: 0.72rem;
    margin-left: 0.3rem; font-family: monospace;
    cursor: help;
}
.risk-classes code {
    background: #fff3cd; padding: 0.15rem 0.4rem;
    border-radius: 3px; margin: 0 0.2rem;
}
.risk-model-panel {
    background: #fafbfc; border: 1px solid #e0e4e8;
    border-radius: 4px; padding: 0.8rem 1rem;
}
.trace-cards {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
    gap: 0.8rem; margin: 0.8rem 0;
}
.trace-card {
    background: #f6f8fa; border: 1px solid #d8dde2; border-radius: 4px;
    padding: 0.7rem 0.9rem;
}
.trace-card-value {
    font-size: 1.4rem; font-weight: 700; color: var(--blue);
    line-height: 1.1;
}
.trace-card-label {
    font-size: 0.78rem; font-weight: 600; color: var(--grey);
    text-transform: uppercase; letter-spacing: 0.04em;
    margin-top: 0.3rem;
}
.trace-card-hint {
    font-size: 0.72rem; color: #888; margin-top: 0.2rem;
}
.prov-runtime { white-space: nowrap; font-size: 0.8rem; color: var(--grey); }
.prov-runtime a { color: var(--blue); text-decoration: none; }
.prov-runtime a:hover { text-decoration: underline; }
.trace-status { color: #555; font-size: 0.78rem; }
.amr-class-legend {
    margin: 0.6rem 0; padding: 0.6rem 0.9rem; background: #f6f8fa;
    border: 1px solid #d8dde2; border-radius: 4px; font-size: 0.85rem;
    line-height: 1.9;
}
.amr-class-chip {
    display: inline-block; padding: 0.1rem 0.5rem;
    border-radius: 11px; color: white; font-size: 0.68rem;
    font-weight: 600; margin: 0 0.15rem 0 0;
    text-transform: capitalize; vertical-align: middle;
    white-space: nowrap;
}
.amr-chip {
    display: inline-block; margin: 0.1rem 0.35rem 0.1rem 0;
    padding: 0.05rem 0; line-height: 1.55;
}
.amr-chip a { text-decoration: none; }
.amr-chip a:hover strong, .amr-chip a:hover code {
    text-decoration: underline;
}
.mag-amr-cell {
    max-width: 480px; line-height: 1.7;
    font-size: 0.85rem;
}
@media print {
    .section { break-inside: avoid; }
    body { background: white; }
    .sample-filter, .search-box { display: none; }
    .header { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
    .data-table th { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
}
</style>"""


# ═══════════════════════════════════════════════════════════════
# JS
# ═══════════════════════════════════════════════════════════════

HTML_JS = """<script>
(function() {
    'use strict';

    // ── Table sorting ──
    function initSortableTables() {
        document.addEventListener('click', function(e) {
            var th = e.target.closest('th');
            if (!th) return;
            var table = th.closest('table.sortable');
            if (!table) return;
            var tbody = table.querySelector('tbody');
            if (!tbody) return;
            var rows = Array.from(tbody.querySelectorAll('tr'));
            var col = Array.from(th.parentNode.children).indexOf(th);
            var asc = !th.classList.contains('sort-asc');
            var allTh = table.querySelectorAll('th');
            for (var i = 0; i < allTh.length; i++) {
                allTh[i].classList.remove('sort-asc', 'sort-desc');
            }
            th.classList.add(asc ? 'sort-asc' : 'sort-desc');
            rows.sort(function(a, b) {
                var cellA = a.children[col];
                var cellB = b.children[col];
                if (!cellA || !cellB) return 0;
                var va = cellA.textContent.trim();
                var vb = cellB.textContent.trim();
                var na = parseFloat(va.replace(/[,%x]/g, '').replace(/[^0-9.\\-]/g, ''));
                var nb = parseFloat(vb.replace(/[,%x]/g, '').replace(/[^0-9.\\-]/g, ''));
                if (!isNaN(na) && !isNaN(nb)) {
                    return asc ? na - nb : nb - na;
                }
                return asc ? va.localeCompare(vb) : vb.localeCompare(va);
            });
            for (var j = 0; j < rows.length; j++) {
                tbody.appendChild(rows[j]);
            }
        });
    }

    // ── Search box ──
    function initSearchFilter() {
        var input = document.getElementById('globalSearch');
        if (!input) return;
        input.addEventListener('input', function() {
            var filter = this.value.toLowerCase();
            var tables = document.querySelectorAll('table.filterable');
            for (var t = 0; t < tables.length; t++) {
                var trs = tables[t].querySelectorAll('tbody tr');
                for (var i = 0; i < trs.length; i++) {
                    var text = trs[i].textContent.toLowerCase();
                    trs[i].style.display =
                        (filter === '' || text.indexOf(filter) !== -1) ? '' : 'none';
                }
            }
        });
    }

    // ── Sample filter checkboxes ──
    window.toggleSample = function(sampleName, checked) {
        var tds = document.querySelectorAll('td[data-sample="' + sampleName + '"]');
        for (var i = 0; i < tds.length; i++) {
            var row = tds[i].closest('tr');
            if (row) {
                row.style.display = checked ? '' : 'none';
            }
        }
    };

    // ── Smooth scroll for TOC ──
    function initTocScroll() {
        var links = document.querySelectorAll('.toc a');
        for (var i = 0; i < links.length; i++) {
            links[i].addEventListener('click', function(e) {
                var href = this.getAttribute('href');
                if (href && href.startsWith('#')) {
                    e.preventDefault();
                    var target = document.getElementById(href.substring(1));
                    if (target) {
                        target.scrollIntoView({behavior: 'smooth', block: 'start'});
                    }
                }
            });
        }
    }

    document.addEventListener('DOMContentLoaded', function() {
        initSortableTables();
        initSearchFilter();
        initTocScroll();
    });
})();
</script>"""


# ═══════════════════════════════════════════════════════════════
# HTML BUILDERS
# ═══════════════════════════════════════════════════════════════

def build_sample_checkboxes(surv):
    if surv.empty:
        return ''
    items = []
    for s in sorted(surv['Sample'].unique()):
        if 'unclassified' in s.lower():
            continue
        items.append(
            f'<label><input type="checkbox" checked '
            f'onchange="toggleSample(\'{s}\', this.checked)"> {s}</label>')
    return ' '.join(items)


def build_header(org_info, logo_b64, title, subtitle, date_str):
    logo_html = (
        f'<img src="{logo_b64}" class="header-logo" alt="Logo">'
        if logo_b64 else '')
    org_line = f'{org_info.get("name", "")} | {org_info.get("department", "")}'
    return f"""<div class="header">
        {logo_html}
        <div class="header-text">
            <h1>{title}</h1>
            <div class="subtitle">{subtitle}</div>
            <div class="org-info">{org_line} | {date_str}</div>
        </div>
    </div>"""


def build_footer(org_info, date_str):
    return f"""<div class="footer">
        EpiTaxMAG Pipeline v2.0 | {org_info.get('name', '')} -
        {org_info.get('department', '')} | {org_info.get('address', '')} |
        {org_info.get('contact', '')} | {date_str}
    </div>"""


def build_toc():
    """Table of contents organised by evidence layer (reads / contigs /
    MAGs / cross-cutting / provenance) rather than by pipeline phase.
    Each layer is rendered as a header in the sidebar so the reader can
    navigate the report by where the evidence comes from."""
    groups = [
        ('Run overview', [
            ('sec-dag',       'Pipeline workflow'),
            ('sec-runstats',  'Run statistics'),
            ('sec-retention', 'Read retention (raw → filtered)'),
            ('sec-summary',   'Per-sample summary'),
        ]),
        ('Reads layer', [
            ('sec-screening',    'Cross-contamination (FastQ Screen + Sylph)'),
            ('sec-diversity',    'Diversity (rarefaction + Whittaker)'),
            ('sec-phenotypic',   'Whole-genome verification (minimap2)'),
            ('sec-amr-reads',    'AMR in reads (KMA / ResFinder)'),
        ]),
        ('Contigs / assembly', [
            ('sec-assembly', 'Assembly and binning'),
            ('sec-plasmids', 'Plasmids and AMR carriage'),
        ]),
        ('MAGs layer', [
            ('sec-mags',         'MAG catalogue (CheckM2 + Bakta + AMR count)'),
            ('sec-taxonomy',     'Per-sample taxonomic composition'),
            ('sec-amr',          'AMR in MAGs (AMRFinderPlus)'),
            ('sec-heatmap',      'AMR heatmap (organism × class)'),
            ('sec-vfdb',         'Virulence factors (VFDB)'),
            ('sec-integrons',    'Integrons (IntegronFinder)'),
        ]),
        ('Cross-cutting', [
            ('sec-concordance', 'Concordance: reads vs MAG contigs'),
            ('sec-risk',        'Risk model and assessment'),
        ]),
        ('Provenance', [
            ('sec-software',  'Software & versions'),
            ('sec-execution', 'Pipeline execution & performance'),
            ('sec-runparams', 'Run parameters (params.json)'),
        ]),
    ]
    blocks = []
    for header, links in groups:
        items = '\n'.join(
            f'<a href="#{sid}">{_esc(label)}</a>' for sid, label in links)
        blocks.append(
            f'<div class="toc-group-header">{_esc(header)}</div>{items}')
    return f'<div class="toc">{"".join(blocks)}</div>'


def build_report(data, org_info, logo_b64, run_name, results_dir=''):
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')
    plotly_js = plotly.offline.get_plotlyjs()

    surv = data['survival']
    checkm2 = data['checkm2']
    integ = data['integration']
    conc = data['concordance']
    amr_mags = data['amr_mags']
    run_stats = data.get('run_stats')
    tax_df = data['taxonomy']
    vfdb_df = data['vfdb']
    plasmids_df = data['plasmids']
    bakta_df = data['bakta']
    integron_summary = data.get('integron_summary', pd.DataFrame())
    quast_df = data.get('quast', pd.DataFrame())
    bins_df = data.get('bins_count', pd.DataFrame())
    tax_info = data.get('tax_info') or detect_taxonomy_tool(results_dir)
    run_params = data.get('run_params') or load_run_params(results_dir)
    trace_df = data.get('trace') if data.get('trace') is not None else load_trace(results_dir)
    fqs_detail = data.get('fqs_detail', pd.DataFrame())

    def _prov(section_key):
        return html_provenance_panel(section_key, run_params, trace_df)

    # ── Metrics ──
    # Two ways to summarise retention across the run:
    #   - read-weighted (sum_clean / sum_raw): the honest "data loss"
    #     of the run. A small bad sample does not move it; a big bad
    #     sample does. This is the one shown on the metric card and in
    #     the TOTAL row of the survival table — they should match.
    #   - per-sample arithmetic mean: useful as a quality indicator
    #     when sample sizes are similar; misleading otherwise.
    # We keep both visible (card = weighted; "Mean per sample" tooltip
    # in the survival table = arithmetic) so the discrepancy is
    # interpretable instead of confusing.
    n_samples = len(surv) if not surv.empty else 0
    total_raw_gb = surv['raw_gb'].sum() if not surv.empty else 0
    total_clean_gb = surv['clean_gb'].sum() if not surv.empty else 0
    if not surv.empty and surv['raw_reads'].sum() > 0:
        mean_ret = float(surv['clean_reads'].sum() / surv['raw_reads'].sum() * 100)
    else:
        mean_ret = 0.0
    mean_ret_per_sample = (surv['pct_reads'].mean()
                           if not surv.empty else 0)
    n_mags = len(checkm2) if not checkm2.empty else 0
    n_hq = (len(checkm2[checkm2['Quality'] == 'HQ'])
            if not checkm2.empty else 0)
    n_amr = len(integ) if not integ.empty else 0
    n_plasmid_amr = (len(integ[integ['Location'] == 'PLASMID'])
                     if not integ.empty else 0)
    n_critical = (len(integ[integ['Risk_level'] == 'CRITICAL'])
                  if not integ.empty else 0)
    n_vf = len(vfdb_df) if not vfdb_df.empty else 0

    pct_unclass = run_stats['pct_unclassified'] if run_stats else 0
    total_run_reads = run_stats['total_reads'] if run_stats else 0
    unclass_reads = run_stats['unclassified_reads'] if run_stats else 0
    unclass_color = (RED if pct_unclass >= 25
                     else (ORANGE if pct_unclass >= 15 else GREEN))

    metrics_html = f"""<div class="metrics-row">
        <div class="metric-card">
            <div class="value">{n_samples}</div>
            <div class="label">Samples</div>
        </div>
        <div class="metric-card">
            <div class="value">{total_raw_gb:.1f}</div>
            <div class="label">Raw Gb</div>
        </div>
        <div class="metric-card">
            <div class="value">{total_clean_gb:.1f}</div>
            <div class="label">Clean Gb</div>
        </div>
        <div class="metric-card" title="Read-weighted retention across the run (sum of clean reads / sum of raw reads). Per-sample arithmetic mean is {mean_ret_per_sample:.1f}%.">
            <div class="value">{mean_ret:.1f}%</div>
            <div class="label">Run retention<br><small style="font-weight:normal;color:#888">read-weighted</small></div>
        </div>
        <div class="metric-card">
            <div class="value" style="color:{unclass_color}">{pct_unclass}%</div>
            <div class="label">Unclassified</div>
        </div>
        <div class="metric-card">
            <div class="value">{n_mags}</div>
            <div class="label">MAGs ({n_hq} HQ)</div>
        </div>
        <div class="metric-card">
            <div class="value">{n_amr}</div>
            <div class="label">AMR Genes</div>
        </div>
        <div class="metric-card">
            <div class="value" style="color:{RED if n_critical > 0 else GREEN}">
                {n_critical}
            </div>
            <div class="label">Critical Risk</div>
        </div>
    </div>"""

    # ── Alerts ──
    alerts = ''
    if n_critical > 0 and not integ.empty:
        alert_items = ''
        for _, r in integ[integ['Risk_level'] == 'CRITICAL'].iterrows():
            alert_items += (
                f'<p><strong>{r.get("Gene", "")}</strong> '
                f'({r.get("AMR_class", "")}) in '
                f'<em>{r.get("Organism", "")}</em> - '
                f'plasmid score {r.get("Plasmid_score", 0):.2f} '
                f'[{r.get("Sample", "")}]</p>')
        alerts = (
            f'<div class="alert-box">'
            f'<p><strong>ALERT: {n_critical} critical AMR genes on '
            f'plasmids (transferable)</strong></p>{alert_items}</div>')

    # ── Charts ──
    retention_chart = fig_to_html(chart_retention(surv))
    summary_chart = fig_to_html(chart_summary_bar(data))
    taxonomy_chart = fig_to_html(chart_taxonomy_stacked(tax_df, checkm2, tax_info))
    heatmap_chart = fig_to_html(chart_amr_heatmap(integ))
    amr_gene_heatmap = fig_to_html(chart_amr_gene_heatmap(integ, amr_mags))
    risk_chart = fig_to_html(chart_risk_summary(integ))
    bins_chart = fig_to_html(chart_bins_per_binner(bins_df))

    # ── Screening charts ──
    fqs_df = load_fastqscreen_raw(results_dir)
    screening_fqs_html = fig_to_html(build_fastqscreen_100pct(fqs_df))
    target_organisms = [t.strip() for t in
                        "Escherichia,Salmonella,Shigella,Vibrio,Campylobacter,Aeromonas,Arcobacter,Aliarcobacter,Pseudomonas,Legionella,Acinetobacter,Mycobacterium,Klebsiella,Enterococcus,Leptospira,Burkholderia,Phocaeicola,Bacteroides,Dorea,Naegleria,Cryptosporidium,Giardia".split(',')]
    sylph_path = f'{results_dir}/10_tax_sylph/sylph_profile_all.tsv'
    sylph_full = load_sylph_profile(sylph_path)
    screening_sylph_html = fig_to_html(build_sylph_dual_heatmap(sylph_full, target_organisms))

    # ── Run stats table ──
    unclass_warn = ''
    if pct_unclass >= 25:
        unclass_warn = (
            f'<span style="color:{RED}; font-weight:700"> '
            f'(ALERT: &ge;25%)</span>')
    elif pct_unclass >= 15:
        unclass_warn = (
            f'<span style="color:{ORANGE}; font-weight:700"> '
            f'(Warning: &ge;15%)</span>')

    run_stats_table = f"""
    <table class="data-table">
        <thead><tr><th>Metric</th><th>Value</th><th>Assessment</th></tr></thead>
        <tbody>
            <tr><td>Total reads sequenced</td>
                <td><strong>{total_run_reads:,}</strong></td><td>-</td></tr>
            <tr><td>Classified reads (with barcode)</td>
                <td>{total_run_reads - unclass_reads:,}</td><td>-</td></tr>
            <tr><td>Unclassified reads (no barcode)</td>
                <td>{unclass_reads:,}</td>
                <td style="color:{unclass_color}; font-weight:700">
                    {pct_unclass}%{unclass_warn}</td></tr>
            <tr><td>Demultiplexed samples</td>
                <td>{n_samples}</td><td>-</td></tr>
            <tr><td>Total raw Gb (samples)</td>
                <td>{total_raw_gb:.2f}</td><td>-</td></tr>
            <tr><td>Total clean Gb (after filtering)</td>
                <td>{total_clean_gb:.2f}</td><td>-</td></tr>
            <tr><td>Recovered MAGs</td>
                <td>{n_mags} ({n_hq} HQ)</td><td>-</td></tr>
            <tr><td>AMR genes detected</td>
                <td>{n_amr} ({n_plasmid_amr} on plasmids)</td><td>-</td></tr>
            <tr><td>Virulence factors (VFDB)</td>
                <td>{n_vf}</td><td>-</td></tr>
            <tr><td>Critical risks</td>
                <td style="color:{RED if n_critical > 0 else GREEN};
                    font-weight:700">{n_critical}</td>
                <td>{'ALERT' if n_critical > 0 else 'OK'}</td></tr>
        </tbody>
    </table>"""

    header = build_header(
        org_info, logo_b64,
        f'EpiTaxMAG Report v2 - {run_name}',
        f'Comprehensive Report | {n_samples} samples | {n_mags} MAGs | '
        f'{n_amr} AMR genes',
        date_str)
    footer = build_footer(org_info, date_str)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EpiTaxMAG Report v2 - {run_name}</title>
    <script type="text/javascript">{plotly_js}</script>
    {HTML_CSS}
</head>
<body>
    {header}
    <div class="container">
        {metrics_html}
        {alerts}

        {build_toc()}

        <div class="sample-filter" id="sampleFilter">
            <strong>Filter samples:</strong>
            {build_sample_checkboxes(surv)}
            <p class="note" style="margin-top:0.4rem">
                Note: sample filters affect only tables,
                not Plotly charts.</p>
        </div>

        <div class="search-box">
            <input type="text" id="globalSearch"
                   placeholder="Search across all tables: gene, organism, AMR class, sample...">
        </div>

        <!-- ═════════════════ RUN OVERVIEW ═════════════════ -->
        <h2 class="layer-header">Run overview</h2>

        <div class="section" id="sec-dag">
            <div class="section-header">Pipeline workflow</div>
            <div class="section-body">
                <p>End-to-end DAG with the tool and version that produces
                each numbered output directory. Use this as a map; every
                downstream section drills into one branch of the DAG.</p>
                {html_pipeline_dag(tax_info)}
            </div>
        </div>

        <div class="section" id="sec-runstats">
            <div class="section-header">Run statistics</div>
            <div class="section-body">
                <p>Global Nanopore run metrics. The <strong>unclassified</strong>
                read percentage (no barcode) is a direct measure of
                demultiplexing efficiency, not of biological diversity.</p>
                {run_stats_table}
                {_prov('qc-survival')}
            </div>
        </div>

        <div class="section" id="sec-retention">
            <div class="section-header">Read retention (raw → filtered)</div>
            <div class="section-body">
                <p>Raw vs clean by sample after Porechop ABI (adapter trimming) +
                Chopper (quality and length filtering). Highlighted: ≥1 Gb clean
                yield. Below-threshold samples are not necessarily failures —
                they may simply be low-biomass.</p>
                {html_survival_table(surv)}
                {retention_chart}
                {_prov('qc-survival')}
            </div>
        </div>

        <div class="section" id="sec-summary">
            <div class="section-header">Per-sample summary</div>
            <div class="section-body">
                <p>Per-sample comparative count of MAGs, AMR genes, virulence
                factors, plasmids and integrons. Use it as the
                index of subsequent sections.</p>
                {summary_chart}
            </div>
        </div>

        <!-- ═════════════════ READS LAYER ═════════════════ -->
        <h2 class="layer-header">Reads layer · evidence directly from the FASTQs</h2>

        <div class="section" id="sec-screening">
            <div class="section-header">Cross-contamination screening (FastQ Screen + Sylph)</div>
            <div class="section-body">
                <p><strong>FastQ Screen</strong> (chart below) maps a subsample
                of reads against a curated genome panel; 100% bars with mutually
                exclusive categories — anything that lights up host or vector
                here means real contamination, not biology.<br>
                <strong>Sylph early-alert</strong> (next chart) shows
                abundance and coverage for the surveillance organism
                list; toggle genus/species and abundance/coverage.</p>
                {screening_fqs_html}
                {screening_sylph_html}
                <p class="note">Surveillance list configurable in
                <code>params.target_organisms</code>.</p>

                <h4 style="margin-top:1.5em">FastQ Screen executive view · one row per sample</h4>
                <p>Single-row-per-sample summary so a 10–20 sample run can
                be scanned at a glance: best surveillance-target hit,
                strongest contamination signal, and a Verdict badge.
                Open the drill-down below to inspect any sample × panel-genome
                cell in detail.</p>
                {html_fastqscreen_summary(fqs_detail, run_params.get('params', {}).get('target_organisms', ''))}
                <details class="prov-panel" style="margin-top:0.6rem">
                    <summary><span class="prov-icon">+</span>
                    <span>Drill-down: full per-sample × panel-genome FastQ Screen table</span>
                    </summary>
                    <p style="margin-top:0.6rem">
                    Five mutually exclusive read-mapping categories taken
                    directly from <code>02_fastqscreen_raw/*_screen.txt</code>;
                    the Category column is an inference of this report
                    (see footnote). Use the global search bar to filter by
                    sample or genome.
                    </p>
                    {html_fastqscreen_detail(fqs_detail, run_params.get('params', {}).get('target_organisms', ''))}
                </details>

                <h4 style="margin-top:1.5em">Bracken vs Sylph at genus level · executive view</h4>
                <p>Per-sample shared / discordant counts and the three biggest
                abundance disagreements between classifiers. <em>B</em> = Bracken-only,
                <em>S</em> = Sylph-only. Open the drill-down for the full
                per-genus table.</p>
                {html_sylph_summary(data['sylph_concordance'])}
                <details class="prov-panel" style="margin-top:0.6rem">
                    <summary><span class="prov-icon">+</span>
                    <span>Drill-down: full Bracken vs Sylph per-genus table</span>
                    </summary>
                    <p style="margin-top:0.6rem">
                    Where Sylph adds (or fails to add) signal over Kraken2/Bracken
                    in this run, gene by gene.
                    </p>
                    {html_sylph_concordance(data['sylph_concordance'])}
                </details>

                {_prov('screening')}
            </div>
        </div>

        <div class="section" id="sec-diversity">
            <div class="section-header">Reads-level diversity (rarefaction + Whittaker)</div>
            <div class="section-body">
                <p>Two complementary views built from the Bracken genus
                counts. <strong>Rarefaction</strong> (left) shows the
                expected number of distinct genera as the read budget
                grows — a curve that has clearly plateaued means we are
                not missing genera by lack of depth; a curve still
                climbing means more sequencing would still discover
                new taxa. The Hurlbert closed form is used (no Monte
                Carlo, deterministic). <strong>Whittaker rank-abundance</strong>
                (right) shows the slope of community dominance — a
                steep curve means a few genera dominate, a shallow one
                means an even community. Compare samples at a glance.</p>
                {fig_to_html(chart_rarefaction(data['rarefaction']))}
                {fig_to_html(chart_whittaker(data['bracken_genus']))}
                {_prov('reads-tax')}
            </div>
        </div>

        <div class="section" id="sec-phenotypic">
            <div class="section-header">Whole-genome verification (minimap2 phenotypic targets)</div>
            <div class="section-body">
                <p>Filtered reads are mapped against the type-strain
                genomes of the surveillance panel
                (<code>params.phenotypic_targets</code>). The triage
                table below interprets <em>mean depth</em> and
                <em>breadth</em> together: a real organism shows depth
                AND breadth aligned, while chance contact gives high
                breadth at 1× but vanishes at 10× and 30×. This is the
                independent confirmation layer for the surveillance
                organisms — orthogonal to Kraken2/Sylph (which classify
                reads against general DBs) and to the MAG layer (which
                needs binnable assemblies).</p>
                {html_phenotypic_table(data.get('phenotypic'))}
                {_prov('phenotypic')}
            </div>
        </div>

        <div class="section" id="sec-amr-reads">
            <div class="section-header">AMR in reads (KMA / ResFinder)</div>
            <div class="section-body">
                <p>KMA hits at the read level against the ResFinder template
                database. High coverage with low depth typically means a
                short-template artifact — depth is the more honest signal
                for read-level resistance.</p>
                {html_kma_reads_table(data['kma'])}
                {_prov('amr-reads')}
            </div>
        </div>

        <!-- ═════════════════ CONTIGS / ASSEMBLY ═════════════════ -->
        <h2 class="layer-header">Contigs / assembly layer</h2>

        <div class="section" id="sec-assembly">
            <div class="section-header">Assembly and binning</div>
            <div class="section-body">
                <p>QUAST metrics on the polished MetaFlye assembly, plus bin
                counts per algorithm (MetaBAT2, MaxBin2, SemiBin2, DAS Tool).
                DAS Tool refines the union into a non-redundant set.</p>
                {html_assembly_table(quast_df)}
                {bins_chart}
                {_prov('assembly')}
                {_prov('binning')}
            </div>
        </div>

        <div class="section" id="sec-plasmids">
            <div class="section-header">Plasmids and AMR carriage</div>
            <div class="section-body">
                <p>Contigs classified as plasmidic by geNomad and the AMR
                genes annotated on them. Score ≥0.9 = confirmed plasmid;
                0.7–0.9 = probable. The AMRFinderPlus pass below the table
                is independent (different model) and useful when geNomad's
                built-in <code>amr_genes</code> column misses something.</p>
                {html_plasmid_table(plasmids_df, integ,
                                    data.get('amr_plasmids'),
                                    amrfinder_catalog=data.get('amrfinder_catalog'))}
                <h4 style="margin-top:1.2em">AMRFinderPlus on plasmid sequences</h4>
                {html_amr_on_plasmids(data.get('amr_plasmids'))}
                <h4 style="margin-top:1.2em">MOB-suite typing · per circular plasmid</h4>
                <p>Replicon family, relaxase (MOB) and predicted mobility for
                plasmids that mob_recon could close. Each row is cross-linked
                to AMR genes detected on the same contig (from AMRFinderPlus
                on MAGs and the AMR-pathogen integration table). A conjugative
                plasmid carrying AMR is the highest dissemination-risk event
                this report can flag.</p>
                {html_mobsuite_detail(data.get('mobsuite'), integ, amr_mags, tax_df)}
                {_prov('plasmids')}
                {_prov('mobsuite')}
            </div>
        </div>

        <!-- ═════════════════ MAGs LAYER ═════════════════ -->
        <h2 class="layer-header">MAGs layer · per-genome evidence</h2>

        <div class="section" id="sec-mags">
            <div class="section-header">MAG catalogue</div>
            <div class="section-body">
                <p>MAGs recovered by MetaBAT2 + MaxBin2 + SemiBin2, refined by
                DAS Tool, evaluated by CheckM2 and classified by
                {_esc(tax_info['short'])}. CDS, tRNA, rRNA, ncRNA and CRISPR
                counts come from Bakta; the AMR genes column shows hits
                per MAG from AMRFinderPlus.</p>
                {html_mag_catalog(checkm2, tax_df, bakta_df, amr_mags)}
                {_prov('mags')}
            </div>
        </div>

        <div class="section" id="sec-taxonomy">
            <div class="section-header">Per-sample taxonomic composition (MAGs)</div>
            <div class="section-body">
                <p>{_esc(_taxonomy_section_caption(tax_info))}</p>
                {taxonomy_chart}
                {_taxonomy_assignment_table(tax_df, tax_info)}
                {_prov('mag-taxonomy')}
            </div>
        </div>

        <div class="section" id="sec-amr">
            <div class="section-header">AMR in MAGs (AMRFinderPlus)</div>
            <div class="section-body">
                <p>AMR genes detected by AMRFinderPlus per MAG.
                Chromosome/plasmid location is inferred from geNomad.
                The <em>Confidence</em> column is a 0–100 composite (identity,
                coverage, method, KMA-MAG concordance, MAG quality); the
                <em>Risk</em> badge follows the deterministic rule set documented
                in the “Risk model” section — hover the rule ID for the trigger.</p>
                {html_amr_detail(integ, amr_mags)}
                {_prov('amr-mags')}
            </div>
        </div>

        <div class="section" id="sec-heatmap">
            <div class="section-header">AMR profile heatmap · gene-level + class summary</div>
            <div class="section-body">
                <p>Two complementary views of the same AMR landscape.
                The <strong>gene-level heatmap</strong> below clusters
                organisms by their resistance fingerprint and groups
                columns by antibiotic class — co-occurrence patterns
                across MAGs are immediately visible.
                <strong>Colour = coverage %</strong> (not identity):
                identity saturates above 95% post-Nanopore polishing
                and the 90% filter compresses every cell to the same
                band, while coverage keeps a real gradient from the
                80% threshold up to 100%, distinguishing full-length
                gene calls from edge-of-contig truncations. Identity
                stays in the per-cell hover. The <strong>class
                summary</strong> further down counts genes per
                organism × class for a bird's-eye view of dissemination.</p>
                {html_amr_class_legend()}
                {amr_gene_heatmap}
                <h4 style="margin-top:1em">Class-level summary</h4>
                <p>Intensity = number of AMR genes per organism × class
                combination. Useful when the gene-level heatmap is dense:
                a column lit up across many organisms is a signal of
                horizontal spread; a single dark cell is an organism-specific
                accumulation.</p>
                {heatmap_chart}
                {_prov('amr-heatmap')}
            </div>
        </div>

        <div class="section" id="sec-vfdb">
            <div class="section-header">Virulence factors (ABRicate / VFDB)</div>
            <div class="section-body">
                <p>Virulence factors detected by ABRicate against VFDB
                ({n_vf} total hits). VFDB hits are mostly informational at
                the metagenomic stage — they say <em>what is encoded</em>, not
                what is being expressed in the sample.</p>
                {html_vfdb_table(vfdb_df)}
                {_prov('vfdb')}
            </div>
        </div>

        <div class="section" id="sec-integrons">
            <div class="section-header">Integrons (IntegronFinder)</div>
            <div class="section-body">
                <p>Integron summary. Only MAGs with at least one integron or
                integrase are shown. <code>intI</code> = integrase-encoding,
                evidence of an active capture/excision platform.</p>
                {html_integron_table(integron_summary, amr_mags, tax_df)}
                {_prov('integrons')}
            </div>
        </div>

        <!-- ═════════════════ CROSS-CUTTING ═════════════════ -->
        <h2 class="layer-header">Cross-cutting evidence and risk</h2>

        <div class="section" id="sec-concordance">
            <div class="section-header">Concordance: reads (KMA) vs MAG contigs (AMRFinderPlus)</div>
            <div class="section-body">
                <p>Cross-validation of AMR signals from reads (KMA/ResFinder)
                against assembled contigs (AMRFinderPlus). Detection in both
                = maximum confidence; reads-only typically means an AMR-bearing
                contig that did not bin into a quality MAG; contigs-only
                means insufficient read depth for KMA on that template.</p>
                {html_concordance_table(conc)}
                {_prov('concordance')}
            </div>
        </div>

        <div class="section" id="sec-risk">
            <div class="section-header">Risk model and assessment</div>
            <div class="section-body">
                <p>The pipeline assigns a risk level per AMR call by walking a
                deterministic ordered rule set. Both the rules and the high-risk
                antibiotic-class set are visible below — every call in the AMR
                detail table carries the matching rule ID, so the decision is
                fully traceable.</p>
                {html_risk_model_panel()}
                <h4 style="margin-top:1.2em">Risk distribution for this run</h4>
                {risk_chart}
                {alerts}
                {_prov('risk')}
            </div>
        </div>

        <!-- ═════════════════ PROVENANCE ═════════════════ -->
        <h2 class="layer-header">Provenance</h2>

        <div class="section" id="sec-software">
            <div class="section-header">Software & versions</div>
            <div class="section-body">
                <p>{len(TOOL_REGISTRY)} tools comprise the EpiTaxMAG pipeline.
                The active genome-taxonomy classifier
                (<strong>{_esc(tax_info['label'])}</strong>) is shown in
                its pipeline-order slot. Per-section parameter detail lives
                inside the “Source &amp; parameters” panel of each section.</p>
                {html_software_table(tax_info)}
            </div>
        </div>

        <div class="section" id="sec-execution">
            <div class="section-header">Pipeline execution &amp; performance</div>
            <div class="section-body">
                <p>Resource and timing profile of every Nextflow task in this
                run, parsed from <code>00_reports/trace.txt</code>. The cards
                below give the run-level totals; the per-layer tables break
                tasks down by stage so bottlenecks (Kraken2/Kaiju RAM,
                MetaFlye walltime, Bakta per-bin scaling) are easy to spot.
                Anomalies are flagged automatically against fixed thresholds
                ({TRACE_ANOMALY_RAM_GB} GB peak RSS, {TRACE_ANOMALY_DURATION_S // 60} min realtime).</p>
                {html_trace_summary_cards(trace_summary(trace_df), run_params)}
                <h4 style="margin-top:1.5em">Anomalies (heavy or slow tasks)</h4>
                {html_trace_anomalies(trace_df)}
                <h4 style="margin-top:1.5em">Top 5 by realtime</h4>
                {html_trace_top_table(trace_df, 'Realtime_s', 5, '')}
                <h4 style="margin-top:1.5em">Top 5 by peak RSS</h4>
                {html_trace_top_table(trace_df, 'Peak_rss_b', 5, '')}
                <h4 style="margin-top:1.5em">Resource distribution</h4>
                <p class="note">Violin plots when a process has ≥3 tasks
                (per-bin Bakta/IntegronFinder, NanoPlot raw+filtered, etc.);
                strip plots otherwise. Log scale on the x-axis for realtime
                and peak RSS — tasks span several orders of magnitude.
                Colour: TAX = blue, MAG = orange, SETUP = green.</p>
                {html_trace_distributions(trace_df)}
                <h4 style="margin-top:1.5em">Tasks grouped by layer</h4>
                {html_trace_by_layer(trace_df)}
            </div>
        </div>

        <div class="section" id="sec-runparams">
            <div class="section-header">Run parameters (params.json)</div>
            <div class="section-body">
                <p>Full snapshot of every <code>params.*</code> value used by
                Nextflow for this run, plus run metadata
                (start/end timestamps, duration, profile, command line). This
                table is read from <code>00_reports/params.json</code>, which
                is written by <code>workflow.onComplete</code>; older runs
                without this file render a notice.</p>
                {html_run_parameters_panel(run_params)}
            </div>
        </div>

    </div>
    {footer}
    {HTML_JS}
</body>
</html>"""
    return html


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive EpiTaxMAG HTML report v2')
    parser.add_argument('--results-dir', required=True,
                        help='Path to results directory')
    parser.add_argument('--run-name', required=True,
                        help='Run name identifier')
    parser.add_argument('--input-dir', default=None,
                        help='Dir with original fastq.gz (for unclassified stats)')
    parser.add_argument('--branding', default='branding/',
                        help='Directory with logo.png and org.info')
    parser.add_argument('--output', required=True,
                        help='Output HTML file path')
    args = parser.parse_args()

    rd = args.results_dir
    print(f'[INFO] Generating EpiTaxMAG v2 report for {args.run_name}')

    org_info, logo_b64 = load_branding(args.branding)
    print(f'[INFO] Branding: {org_info.get("name", "")}')

    # Run stats (unclassified)
    run_stats = None
    if args.input_dir:
        print(f'[INFO] Counting reads from {args.input_dir}...')
        run_stats = load_run_stats(args.input_dir)
        print(f'[INFO] Run: {run_stats["total_reads"]:,} total reads, '
              f'{run_stats["unclassified_reads"]:,} unclassified '
              f'({run_stats["pct_unclassified"]}%)')

    # Load all data
    print('[INFO] Loading QC data...')
    qc_raw = load_qc_csv(f'{rd}/01_qc_raw/qc_raw_metrics.csv')
    qc_filt = load_qc_csv(f'{rd}/06_qc_filtered/qc_filtered_metrics.csv')
    survival = compute_survival(qc_raw, qc_filt)

    print('[INFO] Loading CheckM2...')
    checkm2 = load_checkm2(rd)

    tax_info = detect_taxonomy_tool(rd)
    print(f'[INFO] Loading taxonomy from {tax_info["full"]} '
          f'({tax_info["subdir"]})...')
    taxonomy = load_taxonomy(rd, tax_info)

    print('[INFO] Loading AMRFinderPlus MAGs...')
    amr_mags = load_amr_mags(rd)

    print('[INFO] Loading KMA...')
    kma = load_kma(rd)

    print('[INFO] Loading geNomad plasmids...')
    plasmids = load_plasmids(rd)

    print('[INFO] Loading AMRFinderPlus on plasmid contigs...')
    amr_plasmids = load_amr_plasmids(rd)

    print('[INFO] Loading FastQ Screen detail tables...')
    fqs_detail = load_fastqscreen_detail(rd)
    if not fqs_detail.empty:
        print(f'       fastqscreen detail: {len(fqs_detail)} rows '
              f'across {fqs_detail["Sample"].nunique()} samples')

    print('[INFO] Loading Bracken / Sylph genus tables for concordance...')
    bracken_genus = load_bracken_genus(rd)
    sylph_genus = load_sylph_genus(rd)
    sylph_concordance = compare_bracken_sylph(bracken_genus, sylph_genus)

    print('[INFO] Computing rarefaction curves...')
    rarefaction_df = compute_rarefaction(bracken_genus)

    print('[INFO] Loading phenotypic minimap2 validation...')
    pheno_df = load_phenotypic_validation(rd)
    if not pheno_df.empty:
        print(f'       phenotypic targets: {len(pheno_df)} sample×species rows')
    if not sylph_concordance.empty:
        print(f'       Bracken-Sylph concordance: {len(sylph_concordance)} rows')

    print('[INFO] Loading Nextflow trace.txt and params.json...')
    trace_df = load_trace(rd)
    run_params = load_run_params(rd)

    # AMRFinderPlus reference-gene catalogue: turn NF identifiers
    # (geNomad's amr_genes column, NF033135 etc.) into human gene
    # name + class + subclass + family description.
    amrfinder_db = (run_params.get('params', {}) or {}).get('amrfinder_db')
    amrfinder_catalog = load_amrfinder_catalog(amrfinder_db)
    if amrfinder_catalog:
        n_distinct = len({k.split('.', 1)[0] for k in amrfinder_catalog})
        print(f'       AMRFinderPlus catalogue: {n_distinct} NF identifiers')
    if not trace_df.empty:
        print(f'       trace rows: {len(trace_df)} '
              f'({int((trace_df["status"] == "CACHED").sum())} cached)')
    if run_params:
        print(f'       params.json: {len(run_params.get("params", {}))} parameters')

    print('[INFO] Loading AMR-pathogen integration...')
    integration = load_integration(rd)

    print('[INFO] Loading VFDB (ABRicate)...')
    vfdb = load_vfdb(rd)

    print('[INFO] Loading Bakta summaries...')
    bakta = load_bakta_summaries(rd)

    print('[INFO] Loading IntegronFinder...')
    integrons_raw = load_integrons(rd)
    integron_summary = summarize_integrons(integrons_raw)

    print('[INFO] Loading MOB-suite...')
    mobsuite = load_mobsuite(rd)

    print('[INFO] Loading QUAST assembly metrics...')
    quast = load_quast(rd)

    print('[INFO] Counting bins per binner...')
    bins_count = count_bins(rd)

    print('[INFO] Computing concordance KMA-MAG...')
    concordance = compute_concordance(kma, amr_mags)

    # Count intI annotations for summary
    n_intI = 0
    if not integrons_raw.empty and 'annotation' in integrons_raw.columns:
        n_intI = int((integrons_raw['annotation'] == 'intI').sum())

    print(f'[INFO] Data summary:')
    print(f'       Samples: {len(survival)}')
    print(f'       MAGs: {len(checkm2)}')
    print(f'       Taxonomy: {len(taxonomy)}')
    print(f'       AMR integration: {len(integration)}')
    print(f'       KMA hits: {len(kma)}')
    print(f'       VFDB hits: {len(vfdb)}')
    print(f'       Plasmids: {len(plasmids)}')
    print(f'       Bakta summaries: {len(bakta)}')
    print(f'       Integron MAGs with hits: {len(integron_summary)}')
    print(f'       Integron intI annotations: {n_intI}')
    print(f'       MOB-suite plasmids: {len(mobsuite)}')
    print(f'       QUAST samples: {len(quast)}')
    print(f'       Bins count rows: {len(bins_count)}')
    print(f'       Concordance rows: {len(concordance)}')
    print(f'       AMR-on-plasmid hits: {len(amr_plasmids)}')

    data = {
        'survival': survival,
        'checkm2': checkm2,
        'taxonomy': taxonomy,
        'amr_mags': amr_mags,
        'kma': kma,
        'integration': integration,
        'concordance': concordance,
        'run_stats': run_stats,
        'vfdb': vfdb,
        'plasmids': plasmids,
        'bakta': bakta,
        'integron_summary': integron_summary,
        'integrons_raw': integrons_raw,
        'mobsuite': mobsuite,
        'quast': quast,
        'bins_count': bins_count,
        'tax_info': tax_info,
        'amr_plasmids': amr_plasmids,
        'fqs_detail': fqs_detail,
        'sylph_concordance': sylph_concordance,
        'bracken_genus': bracken_genus,
        'rarefaction': rarefaction_df,
        'phenotypic': pheno_df,
        'trace': trace_df,
        'run_params': run_params,
        'amrfinder_catalog': amrfinder_catalog,
    }

    # Ensure output directory exists. `args.output` may be a bare
    # filename (no directory part) when the script runs inside a
    # Nextflow work dir; fall back to the current directory in that case.
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)

    print('[INFO] Building HTML report...')
    html = build_report(data, org_info, logo_b64, args.run_name, results_dir=args.results_dir)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f'[OK] Report written: {args.output} ({size_mb:.1f} MB)')
    print('[DONE]')


if __name__ == '__main__':
    main()
