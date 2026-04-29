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

SOFTWARE_VERSIONS = [
    ('NanoPlot', '1.46.2', 'QC reads', 'https://github.com/wdecoster/NanoPlot'),
    ('FastQC', '0.12.1', 'QC reads', 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'),
    ('FastQ Screen', '0.16.0', 'Contamination screening', 'https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/'),
    ('Porechop ABI', '0.5.1', 'Adapter trimming', 'https://github.com/bonsai-team/Porechop_ABI'),
    ('Chopper', '0.12.0', 'Quality/length filtering', 'https://github.com/wdecoster/chopper'),
    ('minimap2', '2.30', 'Mapping and host removal', 'https://github.com/lh3/minimap2'),
    ('samtools', '1.21', 'BAM/SAM manipulation', 'https://www.htslib.org/'),
    ('Kraken2', '2.17.1', 'Read taxonomic classification', 'https://ccb.jhu.edu/software/kraken2/'),
    ('Bracken', '3.1', 'Abundance re-estimation', 'https://ccb.jhu.edu/software/bracken/'),
    ('Kaiju', '1.10.1', 'Protein taxonomic classification', 'https://bioinformatics-centre.github.io/kaiju/'),
    ('Sylph', '0.9.0', 'k-mer taxonomic profiling', 'https://github.com/bluenote-1577/sylph'),
    ('KMA', '1.6.8', 'AMR reads (ResFinder)', 'https://bitbucket.org/genomicepidemiology/kma/'),
    ('MultiQC', '1.33', 'QC aggregation', 'https://multiqc.info/'),
    ('MetaFlye', '2.9.6', 'Metagenomic assembly', 'https://github.com/fenderglass/Flye'),
    ('Medaka', '2.2.1', 'Consensus polishing', 'https://github.com/nanoporetech/medaka'),
    ('QUAST', '5.3.0', 'Assembly QC', 'https://quast.sourceforge.net/'),
    ('MetaBAT2', '2.18', 'Binning (coverage)', 'https://bitbucket.org/berkeleylab/metabat/'),
    ('MaxBin2', '2.2.7', 'Binning (EM)', 'https://sourceforge.net/projects/maxbin2/'),
    ('SemiBin2', '2.2.1', 'Binning (deep learning)', 'https://github.com/BigDataBiology/SemiBin'),
    ('DAS Tool', '1.1.7', 'Bin refinement', 'https://github.com/cmks/DAS_Tool'),
    ('CheckM2', '1.1.0', 'MAG QC', 'https://github.com/chklovski/CheckM2'),
    ('GTDB-Tk', '2.7.0 (r232)', 'Genome taxonomy', 'https://ecogenomics.github.io/GTDBTk/'),
    ('AMRFinderPlus', '4.2.7', 'AMR contigs', 'https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/'),
    ('geNomad', '1.12.0', 'Plasmids/viruses', 'https://github.com/apcamargo/genomad'),
    ('ABRicate', '1.0.1', 'Virulence (VFDB)', 'https://github.com/tseemann/abricate'),
    ('MOB-suite', '3.1.9', 'Plasmid typing', 'https://github.com/phac-nml/mob-suite'),
    ('IntegronFinder', '2.0.5', 'Integrons', 'https://github.com/gem-pasteur/Integron_Finder'),
    ('Bakta', '1.12.0', 'Genome annotation', 'https://github.com/oschwengers/bakta'),
]


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


def load_taxonomy(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/24_taxonomy_gtdbtk/*/*_taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                continue
            for _, r in df.iterrows():
                cls = str(r.iloc[1]) if len(r) > 1 else ''
                levels = {}
                for part in cls.split(';'):
                    p = part.strip()
                    for px, lv in [('d__', 'domain'), ('p__', 'phylum'), ('c__', 'class'),
                                   ('o__', 'order'), ('f__', 'family'), ('g__', 'genus'),
                                   ('s__', 'species')]:
                        if p.startswith(px):
                            levels[lv] = p[len(px):] or ''
                org = f"{levels.get('genus', '')} {levels.get('species', '')}".strip()
                if not org:
                    org = 'Unclassified'
                rows.append({
                    'Sample': sample,
                    'Bin': str(r.iloc[0]).strip(),
                    'Taxonomy': cls,
                    'Organism': org,
                    'Domain': levels.get('domain', ''),
                    'Phylum': levels.get('phylum', ''),
                    'Class': levels.get('class', ''),
                    'Order': levels.get('order', ''),
                    'Family': levels.get('family', ''),
                    'Genus': levels.get('genus', ''),
                    'Species': levels.get('species', ''),
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
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/26_genomad/*/*_plasmid_summary.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                contig = str(r.iloc[0]).strip().split('|')[0]
                score = float(r.get('plasmid_score', r.get('score', 0)))
                n_genes = int(r.get('n_genes', 0))
                amr_genes = str(r.get('amr_genes', ''))
                rows.append({
                    'Sample': sample, 'Contig': contig,
                    'Plasmid_score': score,
                    'Length': int(r.get('length', 0)),
                    'N_genes': n_genes,
                    'AMR_genes': amr_genes,
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


def chart_taxonomy_stacked(tax_df, checkm2_df):
    """Stacked bar: organisms per sample colored by Class (top 10 + Other)."""
    if tax_df.empty:
        return None
    # Merge with checkm2 to get only quality-passed MAGs
    df = tax_df.copy()
    if not checkm2_df.empty:
        valid_bins = set(zip(checkm2_df['Sample'], checkm2_df['Bin']))
        df = df[df.apply(lambda r: (r['Sample'], r['Bin']) in valid_bins, axis=1)]
    if df.empty:
        return None

    # Use Class level for grouping
    df['Color_group'] = df['Class'].apply(
        lambda v: v if v and str(v).strip() else 'Unknown')

    # Keep top 10 classes by count, rest becomes "Other"
    class_counts = df['Color_group'].value_counts()
    top_classes = class_counts.nlargest(10).index.tolist()
    df['Color_group'] = df['Color_group'].apply(
        lambda v: v if v in top_classes else 'Other')

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
        title=dict(text='Per-Sample Taxonomic Composition (Class)',
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

def html_pipeline_dag():
    return """<div class="dag-container"><div class="dag-flow">
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
        <div class="dag-box tax">26 GTDB-Tk 2.7.0 (r232)</div>
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
    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Raw Reads</th><th>Raw Gb</th>'
        '<th>Clean Reads</th><th>Clean Gb</th><th>% Retention</th>'
        f'</tr></thead><tbody>{rows}</tbody></table>')


def html_mag_catalog(checkm2_df, tax_df, bakta_df):
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
    # Merge Bakta summary for tRNA/rRNA
    if not bakta_df.empty:
        bakta_merge = bakta_df[['Sample', 'Bin', 'tRNAs', 'rRNAs']].copy()
        df = df.merge(bakta_merge, on=['Sample', 'Bin'], how='left')
    else:
        df['tRNAs'] = np.nan
        df['rRNAs'] = np.nan
    df = df.fillna('')

    rows = ''
    for _, r in df.iterrows():
        trna_val = str(int(r['tRNAs'])) if r['tRNAs'] != '' and r['tRNAs'] != 0 else '-'
        rrna_val = str(int(r['rRNAs'])) if r['rRNAs'] != '' and r['rRNAs'] != 0 else '-'
        try:
            trna_val = str(int(float(r['tRNAs']))) if r['tRNAs'] != '' else '-'
        except (ValueError, TypeError):
            trna_val = '-'
        try:
            rrna_val = str(int(float(r['rRNAs']))) if r['rRNAs'] != '' else '-'
        except (ValueError, TypeError):
            rrna_val = '-'
        rows += (
            f'<tr><td class="sample" data-sample="{r.get("Sample", "")}">'
            f'{r.get("Sample", "")}</td>'
            f'<td>{_esc(r["Bin"])}</td>'
            f'<td><em>{_esc(r.get("Organism", ""))}</em></td>'
            f'<td>{r["Completeness"]:.1f}%</td>'
            f'<td>{r["Contamination"]:.1f}%</td>'
            f'<td>{QC_BADGE.get(r["Quality"], "")}</td>'
            f'<td>{trna_val}</td>'
            f'<td>{rrna_val}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Bin</th><th>Organism</th>'
        '<th>Completeness</th><th>Contamination</th><th>Quality</th>'
        '<th>tRNAs</th><th>rRNAs</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>')


def html_amr_detail(integ_df, amr_mags_df):
    if integ_df.empty:
        return '<p class="no-data">No AMR genes detected in MAGs.</p>'
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
                f'title="NCBI Protein">{_esc(accession)}</a>')
        if gene_sym:
            card_url = (f'https://card.mcmaster.ca/ontology/query?query='
                        f'{gene_sym}')
            links_parts.append(
                f'<a href="{card_url}" target="_blank" rel="noopener" '
                f'title="CARD">CARD</a>')
        links_html = ' | '.join(links_parts) if links_parts else '-'

        method_val = str(r.get('Method', ''))

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
            f'<td>{RISK_BADGE.get(r.get("Risk_level", ""), "")}</td>'
            f'<td class="links-cell">{links_html}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>Organism</th><th>Gene</th><th>Description</th>'
        '<th>AMR Class</th><th>Location</th>'
        '<th>% Identity</th><th>% Coverage</th><th>Method</th>'
        '<th>Confidence</th><th>Risk</th><th>Links</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Confidence: composite score (0-100) based on identity, '
        'coverage, detection method, KMA-MAG concordance and MAG quality. '
        'Links: NCBI Protein (accession) and CARD (resistance database).</p>')


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


def html_plasmid_table(plasmids_df, integ_df):
    if plasmids_df.empty:
        return '<p class="no-data">No plasmids detected by geNomad.</p>'

    # Also show AMR genes in plasmids
    plasmid_amr_rows = ''
    if not integ_df.empty:
        plas_amr = integ_df[integ_df['Location'] == 'PLASMID'].copy()
    else:
        plas_amr = pd.DataFrame()

    rows = ''
    for s in sorted(plasmids_df['Sample'].unique()):
        s_df = plasmids_df[plasmids_df['Sample'] == s]
        n_plas = len(s_df)
        mean_score = s_df['Plasmid_score'].mean()
        amr_in_plas = []
        if not plas_amr.empty:
            s_amr = plas_amr[plas_amr['Sample'] == s]
            for _, ar in s_amr.iterrows():
                amr_in_plas.append(str(ar.get('Gene', '')))
        amr_str = ', '.join(amr_in_plas) if amr_in_plas else '-'
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{s}">{s}</td>'
            f'<td>{n_plas}</td>'
            f'<td>{mean_score:.3f}</td>'
            f'<td class="desc-cell">{_esc(amr_str)}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>N Plasmids</th><th>Mean Score</th>'
        '<th>AMR Genes on Plasmids</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Plasmid score (geNomad): &ge;0.9 confirmed, '
        '0.7-0.9 probable, &lt;0.7 uncertain. AMR genes on conjugative plasmids '
        'pose a risk of horizontal dissemination.</p>')


def html_integron_table(integron_summary_df):
    if integron_summary_df.empty:
        return '<p class="no-data">No integrons detected.</p>'
    rows = ''
    for _, r in integron_summary_df.iterrows():
        has_intI_badge = (
            f'<span class="badge" style="background:{RED}">Si</span>'
            if r['has_intI'] == 'Yes'
            else f'<span class="badge" style="background:{GREEN}">No</span>')
        rows += (
            f'<tr>'
            f'<td class="sample" data-sample="{r["Sample"]}">{r["Sample"]}</td>'
            f'<td>{_esc(r["Bin"])}</td>'
            f'<td>{r["N_complete"]}</td>'
            f'<td>{r["N_In0"]}</td>'
            f'<td>{r["N_CALIN"]}</td>'
            f'<td>{has_intI_badge}</td></tr>\n')

    return (
        '<table class="data-table sortable filterable"><thead><tr>'
        '<th>Sample</th><th>MAG</th><th>N Complete</th>'
        '<th>N In0</th><th>N CALIN</th><th>Integrase (intI)</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>'
        '<p class="note">Complete integrons contain integrase + cassettes. '
        'In0 = integrase only without cassettes. CALIN = cassettes without integrase. '
        'Presence of intI indicates capacity to acquire new genes.</p>')


def html_software_table():
    rows = ''
    for name, version, desc, url in SOFTWARE_VERSIONS:
        rows += (
            f'<tr><td><strong>{_esc(name)}</strong></td>'
            f'<td>{_esc(version)}</td>'
            f'<td>{_esc(desc)}</td>'
            f'<td><a href="{url}" target="_blank">{_esc(url)}</a></td></tr>\n')
    return (
        '<table class="data-table sortable"><thead><tr>'
        '<th>Software</th><th>Version</th><th>Function</th><th>URL</th>'
        '</tr></thead>'
        f'<tbody>{rows}</tbody></table>')


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
    display: block; padding: 0.3rem 0;
    font-size: 0.9rem;
}
.toc a:hover { text-decoration: underline; }
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
    """Table of contents with anchor links."""
    sections = [
        ('sec-dag', '1. EpiTaxMAG Pipeline - Workflow'),
        ('sec-runstats', '2. Run Statistics'),
        ('sec-retention', '3. Read Retention'),
        ('sec-screening', '4. Screening: FastQ Screen + Early Alert'),
        ('sec-summary', '5. Per-Sample Summary'),
        ('sec-assembly', '6. Assembly and Binning'),
        ('sec-mags', '7. MAG Catalog'),
        ('sec-taxonomy', '8. Taxonomic Composition'),
        ('sec-amr', '9. Antimicrobial Resistance (AMR)'),
        ('sec-heatmap', '10. AMR Heatmap'),
        ('sec-concordance', '11. KMA-MAG Concordance'),
        ('sec-vfdb', '12. Virulence Factors (VFDB)'),
        ('sec-plasmids', '13. Plasmids'),
        ('sec-integrons', '14. Integrons'),
        ('sec-risk', '15. Risk Assessment'),
        ('sec-software', '16. Software Versions'),
    ]
    links = '\n'.join(
        f'<a href="#{sid}">{label}</a>' for sid, label in sections)
    return f'<div class="toc">{links}</div>'


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

    # ── Metrics ──
    n_samples = len(surv) if not surv.empty else 0
    total_raw_gb = surv['raw_gb'].sum() if not surv.empty else 0
    total_clean_gb = surv['clean_gb'].sum() if not surv.empty else 0
    mean_ret = surv['pct_reads'].mean() if not surv.empty else 0
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
        <div class="metric-card">
            <div class="value">{mean_ret:.1f}%</div>
            <div class="label">Mean Retention</div>
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
    taxonomy_chart = fig_to_html(chart_taxonomy_stacked(tax_df, checkm2))
    heatmap_chart = fig_to_html(chart_amr_heatmap(integ))
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

        <!-- 1. Pipeline DAG -->
        <div class="section" id="sec-dag">
            <div class="section-header">
                1. EpiTaxMAG Pipeline - Workflow (31 steps)
            </div>
            <div class="section-body">
                <p>Complete pipeline diagram with software versions.</p>
                {html_pipeline_dag()}
            </div>
        </div>

        <!-- 2. Run Statistics -->
        <div class="section" id="sec-runstats">
            <div class="section-header">
                2. Sequencing Run Statistics
            </div>
            <div class="section-body">
                <p>Global metrics of the Nanopore run.
                The <strong>unclassified</strong> read percentage
                (no barcode) reflects demultiplexing efficiency.</p>
                {run_stats_table}
            </div>
        </div>

        <!-- 3. Read Retention -->
        <div class="section" id="sec-retention">
            <div class="section-header">
                3. Read Retention by Sample
            </div>
            <div class="section-body">
                <p>Raw vs clean comparison after Porechop ABI (trimming) +
                Chopper (Q/length filtering).
                Green cells: &ge;1 Gb clean.</p>
                {html_survival_table(surv)}
                {retention_chart}
            </div>
        </div>

        <!-- 3b. Screening -->
        <div class="section" id="sec-screening">
            <div class="section-header">
                4. Screening — FastQ Screen and Early Alert for Target Organisms
            </div>
            <div class="section-body">
                <p><strong>FastQ Screen</strong> (lower panel): distribution of raw reads mapped
                against reference genomes. 100% bars with mutually exclusive categories.
                <strong>Early Alert</strong> (upper panel): abundance and coverage of surveillance
                panel organisms detected by Sylph. Buttons to toggle genus/species and
                abundance/coverage.</p>
                {screening_fqs_html}
                {screening_sylph_html}
                <p class="note">Organism list configurable in nextflow.config (params.target_organisms).</p>
            </div>
        </div>

        <!-- 4. Summary Bar Plot -->
        <div class="section" id="sec-summary">
            <div class="section-header">
                5. Per-Sample Summary
            </div>
            <div class="section-body">
                <p>Per-sample comparative count: recovered MAGs,
                AMR genes, virulence factors, plasmids and integrons.</p>
                {summary_chart}
            </div>
        </div>

        <!-- 5. Assembly & Binning -->
        <div class="section" id="sec-assembly">
            <div class="section-header">
                6. Assembly and Binning Metrics
            </div>
            <div class="section-body">
                <p>Assembly metrics from QUAST and bin counts per
                binning algorithm (MetaBAT2, MaxBin2, SemiBin2, DAS Tool).</p>
                {html_assembly_table(quast_df)}
                {bins_chart}
            </div>
        </div>

        <!-- 6. MAG Catalog -->
        <div class="section" id="sec-mags">
            <div class="section-header">
                7. MAG Catalog
            </div>
            <div class="section-body">
                <p>MAGs recovered by MetaBAT2 + MaxBin2 + SemiBin2,
                refined by DAS Tool, evaluated by CheckM2 and classified
                by GTDB-Tk (r232). tRNA/rRNA from Bakta.</p>
                {html_mag_catalog(checkm2, tax_df, bakta_df)}
            </div>
        </div>

        <!-- 7. Taxonomy Stacked Bar -->
        <div class="section" id="sec-taxonomy">
            <div class="section-header">
                8. Per-Sample Taxonomic Composition
            </div>
            <div class="section-body">
                <p>Per-sample organism distribution colored by taxonomic
                class (GTDB-Tk r232). Top 10 classes + Other.</p>
                {taxonomy_chart}
            </div>
        </div>

        <!-- 8. AMR Detail -->
        <div class="section" id="sec-amr">
            <div class="section-header">
                9. Antimicrobial Resistance (AMR) in MAGs
            </div>
            <div class="section-body">
                <p>AMR genes detected by AMRFinderPlus. Chromosome/plasmid
                location by geNomad. Composite confidence score
                (identity + coverage + method + concordance + MAG quality).</p>
                {html_amr_detail(integ, amr_mags)}
            </div>
        </div>

        <!-- 9. AMR Heatmap -->
        <div class="section" id="sec-heatmap">
            <div class="section-header">
                10. AMR Heatmap: Organisms vs Resistance Classes
            </div>
            <div class="section-body">
                <p>Intensity = number of AMR genes per organism-class
                combination.</p>
                {heatmap_chart}
            </div>
        </div>

        <!-- 10. Concordance KMA-MAG -->
        <div class="section" id="sec-concordance">
            <div class="section-header">
                11. Concordance: KMA Reads vs MAG Contigs
            </div>
            <div class="section-body">
                <p>Comparison of AMR genes detected in reads (KMA/ResFinder)
                vs assembled contigs (AMRFinderPlus). Detection in both =
                maximum confidence.</p>
                {html_concordance_table(conc)}
            </div>
        </div>

        <!-- 11. Virulence VFDB -->
        <div class="section" id="sec-vfdb">
            <div class="section-header">
                12. Virulence Factors (VFDB)
            </div>
            <div class="section-body">
                <p>Virulence factors detected by ABRicate against the
                VFDB database (Virulence Factor Database).
                {n_vf} total hits.</p>
                {html_vfdb_table(vfdb_df)}
            </div>
        </div>

        <!-- 12. Plasmids -->
        <div class="section" id="sec-plasmids">
            <div class="section-header">
                13. Plasmids (geNomad)
            </div>
            <div class="section-body">
                <p>Contigs classified as plasmidic by geNomad, with
                associated AMR genes. Score &ge;0.9 = confirmed plasmid.</p>
                {html_plasmid_table(plasmids_df, integ)}
            </div>
        </div>

        <!-- 13. Integrons -->
        <div class="section" id="sec-integrons">
            <div class="section-header">
                14. Integrons (IntegronFinder)
            </div>
            <div class="section-body">
                <p>Summary of integrons detected by IntegronFinder.
                Only MAGs with at least one integron or integrase are shown.</p>
                {html_integron_table(integron_summary)}
            </div>
        </div>

        <!-- 14. Risk Assessment -->
        <div class="section" id="sec-risk">
            <div class="section-header">
                15. AMR Risk Assessment
            </div>
            <div class="section-body">
                {risk_chart}
                <div class="note" style="margin-top:1rem">
                    <strong>Risk levels:</strong><br>
                    CRITICAL = key resistance (beta-lactam, quinolone,
                    colistin) on plasmid<br>
                    HIGH = any AMR on plasmid (transferable)<br>
                    MEDIUM = key resistance on chromosome<br>
                    LOW = non-critical AMR on chromosome
                </div>
                {alerts}
            </div>
        </div>

        <!-- 15. Software Versions -->
        <div class="section" id="sec-software">
            <div class="section-header">
                16. Software Versions
            </div>
            <div class="section-body">
                <p>{len(SOFTWARE_VERSIONS)} tools used in the
                EpiTaxMAG pipeline.</p>
                {html_software_table()}
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

    print('[INFO] Loading GTDB-Tk taxonomy...')
    taxonomy = load_taxonomy(rd)

    print('[INFO] Loading AMRFinderPlus MAGs...')
    amr_mags = load_amr_mags(rd)

    print('[INFO] Loading KMA...')
    kma = load_kma(rd)

    print('[INFO] Loading geNomad plasmids...')
    plasmids = load_plasmids(rd)

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
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print('[INFO] Building HTML report...')
    html = build_report(data, org_info, logo_b64, args.run_name, results_dir=args.results_dir)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f'[OK] Report written: {args.output} ({size_mb:.1f} MB)')
    print('[DONE]')


if __name__ == '__main__':
    main()
