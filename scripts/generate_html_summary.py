#!/usr/bin/env python3
"""
generate_html_summary.py
─────────────────────────────────────────────────────────────────────
Genera un reporte HTML interactivo con resumen completo del pipeline
metagenómico Nanopore: métricas de supervivencia, taxonomía integrada,
viabilidad de MAGs, curvas de rarefacción y conclusión logística.

Uso:
    python3 generate_html_summary.py \
        --qc-raw qc_raw_metrics.csv \
        --qc-filtered qc_filtered_metrics.csv \
        --sylph sylph_profile_all.tsv \
        --run-name EPIM232 \
        --output final_report.html

FISABIO — EPIMOL
"""

import argparse
import os
import glob
import sys
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
from datetime import datetime

warnings.filterwarnings('ignore')

# ─── Paleta FISABIO ─────────────────────────────────────────────────

PALETTE = [
    '#2C5F8A', '#E05C2A', '#28A745', '#9B59B6', '#F39C12',
    '#1ABC9C', '#E74C3C', '#3498DB', '#95A5A6', '#D35400',
    '#2ECC71', '#8E44AD', '#16A085', '#C0392B', '#5B9ABF',
]
BLUE   = '#2C5F8A'
ORANGE = '#E05C2A'
GREEN  = '#28A745'
RED    = '#E74C3C'
GREY   = '#6C757D'


# ════════════════════════════════════════════════════════════════════
# PARSERS
# ════════════════════════════════════════════════════════════════════

def parse_qc_csv(path):
    """Parse QC CSV con autodetección de separador y decimal europeo."""
    with open(path) as f:
        header = f.readline()
    sep = '\t' if '\t' in header else (';' if ';' in header else ',')
    df = pd.read_csv(path, sep=sep)
    for col in df.columns:
        if col == 'Sample':
            continue
        # Convertir todo a string, reemplazar coma decimal, y forzar numérico
        df[col] = pd.to_numeric(
            df[col].astype(str).str.replace(',', '.'), errors='coerce'
        )
    return df


def parse_kraken2_report(path):
    """Parse Kraken2 standard report (6 columnas TSV)."""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 6:
                try:
                    rows.append({
                        'pct': float(parts[0].strip()),
                        'reads_clade': int(parts[1].strip()),
                        'reads_direct': int(parts[2].strip()),
                        'rank': parts[3].strip(),
                        'taxid': parts[4].strip(),
                        'name': parts[5].strip()
                    })
                except ValueError:
                    continue
    return pd.DataFrame(rows)


def parse_bracken_species(path):
    """Parse Bracken species-level output (.bracken.S.txt)."""
    try:
        return pd.read_csv(path, sep='\t')
    except Exception:
        return pd.DataFrame()


def parse_kaiju_species(path):
    """Parse Kaiju species table (kaiju2table output)."""
    try:
        df = pd.read_csv(path, sep='\t')
        df.columns = [c.strip() for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()


def parse_sylph_profile(path):
    """Parse Sylph profile TSV y clasifica por DB."""
    df = pd.read_csv(path, sep='\t')

    def extract_species(cn):
        parts = str(cn).split()
        if len(parts) < 3:
            return str(cn)[:40]
        skip = 1
        if len(parts) > skip and parts[skip].upper() == 'MAG:':
            skip += 1
        genus = parts[skip] if skip < len(parts) else ''
        species = parts[skip + 1] if skip + 1 < len(parts) else ''
        return f'{genus} {species}'.strip() or str(cn)[:40]

    def extract_sample(s):
        base = os.path.basename(str(s))
        for ext in ['.clean.fastq.gz', '.filtered.fastq.gz', '.fastq.gz', '.fastq', '.fq.gz']:
            base = base.replace(ext, '')
        return base

    def classify_db(gf):
        gl = str(gf).lower()
        if 'imgvr' in gl:
            return 'Virus'
        elif 'fungi' in gl:
            return 'Fungi'
        return 'Bacteria/Archaea'

    df['Species'] = df['Contig_name'].apply(extract_species)
    df['Sample'] = df['Sample_file'].apply(extract_sample)
    df['DB'] = df['Genome_file'].apply(classify_db)
    df['Tax_abund'] = pd.to_numeric(df['Taxonomic_abundance'], errors='coerce').fillna(0)
    df['Eff_cov'] = pd.to_numeric(df['Eff_cov'], errors='coerce').fillna(0)
    return df


def parse_fastqscreen_results(screen_files):
    """Parse FastQ Screen _screen.txt files (todos los del directorio de trabajo)."""
    all_data = []
    for fpath in sorted(screen_files):
        sample = os.path.basename(fpath).replace('_screen.txt', '')
        try:
            with open(fpath) as f:
                lines = f.readlines()
        except Exception:
            continue
        header_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith('Genome'):
                header_idx = i
                break
        if header_idx is None:
            continue
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if not line or line.startswith('%') or line.startswith('#'):
                break
            parts = line.split('\t')
            if len(parts) >= 12:
                genome = parts[0]
                try:
                    reads_processed = int(parts[1])
                    pct_one_hit_one = float(parts[5])
                    pct_multi_hit_one = float(parts[7])
                    pct_one_hit_multi = float(parts[9])
                    pct_multi_hit_multi = float(parts[11])
                except (ValueError, IndexError):
                    continue
                total_mapped = pct_one_hit_one + pct_multi_hit_one + pct_one_hit_multi + pct_multi_hit_multi
                all_data.append({
                    'Sample': sample,
                    'Genome': genome,
                    'Reads_processed': reads_processed,
                    'Pct_mapped': round(total_mapped, 2),
                    'Pct_unique': round(pct_one_hit_one, 2),
                    'Pct_multi': round(pct_multi_hit_one + pct_one_hit_multi + pct_multi_hit_multi, 2)
                })
    return pd.DataFrame(all_data) if all_data else pd.DataFrame()


def parse_kma_results(res_files):
    """Parse KMA .res result files."""
    all_data = []
    for fpath in sorted(res_files):
        sample = os.path.basename(fpath).replace('.res', '')
        try:
            with open(fpath) as f:
                header = f.readline()
            if not header.strip():
                continue
            df = pd.read_csv(fpath, sep='\t')
            if df.empty:
                continue
            cols = df.columns.tolist()
            for _, row in df.iterrows():
                template = str(row.iloc[0]).strip()
                if not template or template.startswith('#'):
                    continue
                identity = 0
                coverage = 0
                depth = 0
                score = 0
                for c in cols:
                    cl = c.strip().lower()
                    if 'template_identity' in cl:
                        identity = float(row[c])
                    elif 'template_coverage' in cl:
                        coverage = float(row[c])
                    elif cl == 'depth':
                        depth = float(row[c])
                    elif cl == 'score':
                        score = float(row[c])
                all_data.append({
                    'Sample': sample,
                    'Gene': template,
                    'Template_Identity': identity,
                    'Template_Coverage': coverage,
                    'Depth': depth,
                    'Score': score
                })
        except Exception:
            continue
    return pd.DataFrame(all_data) if all_data else pd.DataFrame()


def parse_phenotypic_results(path):
    """Parse combined phenotypic coverage TSV from minimap2 verification."""
    if not path or not os.path.isfile(path):
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t')


# ════════════════════════════════════════════════════════════════════
# ANÁLISIS
# ════════════════════════════════════════════════════════════════════

def compute_survival(qc_raw, qc_filtered):
    """Tabla comparativa: reads y Gb antes/después de filtrado."""
    raw = qc_raw[['Sample', 'total_reads', 'mean_length']].copy()
    raw.columns = ['Sample', 'raw_reads', 'raw_mean_len']
    raw['raw_gb'] = raw['raw_reads'] * raw['raw_mean_len'] / 1e9

    filt = qc_filtered[['Sample', 'total_reads', 'mean_length']].copy()
    filt.columns = ['Sample', 'clean_reads', 'clean_mean_len']
    filt['clean_gb'] = filt['clean_reads'] * filt['clean_mean_len'] / 1e9

    df = pd.merge(raw, filt, on='Sample', how='outer').fillna(0)
    df['pct_reads'] = np.where(
        df['raw_reads'] > 0,
        df['clean_reads'] / df['raw_reads'] * 100, 0
    ).round(1)
    df['pct_gb'] = np.where(
        df['raw_gb'] > 0,
        df['clean_gb'] / df['raw_gb'] * 100, 0
    ).round(1)
    return df.sort_values('Sample').reset_index(drop=True)


def compute_mag_viability(sylph_df, survival_df, genome_size_mb, min_coverage):
    """Estima cobertura por organismo y evalúa viabilidad para MAGs."""
    results = []
    bact = sylph_df[sylph_df['DB'] == 'Bacteria/Archaea'].copy()

    for sample in survival_df['Sample'].unique():
        row_surv = survival_df[survival_df['Sample'] == sample]
        if row_surv.empty:
            continue
        clean_gb = row_surv['clean_gb'].values[0]

        s_bact = bact[bact['Sample'] == sample]
        if s_bact.empty:
            continue

        grp = s_bact.groupby('Species').agg(
            abund=('Tax_abund', 'sum'),
            eff_cov=('Eff_cov', 'mean')
        ).sort_values('abund', ascending=False).reset_index()

        for _, row in grp.iterrows():
            abund_frac = row['abund'] / 100.0
            # clean_Mb * fracción_abundancia / genoma_Mb = cobertura estimada
            est_cov = clean_gb * 1000 * abund_frac / genome_size_mb

            results.append({
                'Sample': sample,
                'Species': row['Species'],
                'Abund_pct': round(row['abund'], 2),
                'Eff_cov_sylph': round(row['eff_cov'], 1),
                'Est_cov': round(est_cov, 1),
                'Clean_Gb': round(clean_gb, 3),
                'Viable': row['eff_cov'] >= min_coverage
            })

    return pd.DataFrame(results) if results else pd.DataFrame()


def simulate_rarefaction(abundances, n_steps=50, n_rep=10, seed=42):
    """Simula curva de rarefacción desde abundancias relativas."""
    rng = np.random.default_rng(seed)
    props = abundances / abundances.sum()
    n_total = len(props)
    n_max = 100_000
    steps = np.unique(np.logspace(1, np.log10(n_max), n_steps).astype(int))

    means, lows, highs = [], [], []
    for n in steps:
        counts = [np.sum(rng.multinomial(n, props) > 0) for _ in range(n_rep)]
        means.append(np.mean(counts))
        lows.append(np.percentile(counts, 10))
        highs.append(np.percentile(counts, 90))

    return steps, np.array(means), np.array(lows), np.array(highs), n_total


# ════════════════════════════════════════════════════════════════════
# GRÁFICOS PLOTLY
# ════════════════════════════════════════════════════════════════════

def chart_survival(survival_df):
    """Gráfico comparativo raw vs clean (reads y Gb)."""
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=['Reads', 'Gigabases (Gb)'],
        horizontal_spacing=0.12
    )
    samples = survival_df['Sample'].values

    fig.add_trace(go.Bar(
        name='Raw', x=samples, y=survival_df['raw_reads'],
        marker_color=BLUE, opacity=0.7,
        text=[f'{v:,.0f}' for v in survival_df['raw_reads']],
        textposition='outside', textfont_size=9
    ), row=1, col=1)
    fig.add_trace(go.Bar(
        name='Clean', x=samples, y=survival_df['clean_reads'],
        marker_color=GREEN, opacity=0.85,
        text=[f'{v:,.0f}' for v in survival_df['clean_reads']],
        textposition='outside', textfont_size=9
    ), row=1, col=1)

    fig.add_trace(go.Bar(
        name='Raw Gb', x=samples, y=survival_df['raw_gb'],
        marker_color=BLUE, opacity=0.7, showlegend=False,
        text=[f'{v:.2f}' for v in survival_df['raw_gb']],
        textposition='outside', textfont_size=9
    ), row=1, col=2)
    fig.add_trace(go.Bar(
        name='Clean Gb', x=samples, y=survival_df['clean_gb'],
        marker_color=GREEN, opacity=0.85, showlegend=False,
        text=[f'{v:.2f}' for v in survival_df['clean_gb']],
        textposition='outside', textfont_size=9
    ), row=1, col=2)

    fig.update_layout(
        barmode='group', height=420,
        margin=dict(l=60, r=30, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.04, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif'),
    )
    fig.update_xaxes(tickangle=-40, tickfont_size=9)
    return fig


def chart_taxonomy_sample(sample, bracken_df=None, kraken_df=None, kaiju_df=None,
                          top_n=10):
    """Barras horizontales: top especies de cada herramienta para una muestra."""
    data = []

    if bracken_df is not None and not bracken_df.empty:
        if 'fraction_total_reads' in bracken_df.columns:
            top = bracken_df.nlargest(top_n, 'fraction_total_reads')
            for _, row in top.iterrows():
                data.append({
                    'Species': str(row['name'])[:40],
                    'Tool': 'Bracken',
                    'Abundance': row['fraction_total_reads'] * 100
                })

    if kraken_df is not None and not kraken_df.empty:
        species = kraken_df[kraken_df['rank'] == 'S'].nlargest(top_n, 'pct')
        for _, row in species.iterrows():
            if row['name'].strip().lower() == 'unclassified':
                continue
            data.append({
                'Species': row['name'].strip()[:40],
                'Tool': 'Kraken2',
                'Abundance': row['pct']
            })

    if kaiju_df is not None and not kaiju_df.empty:
        pct_col = 'percent' if 'percent' in kaiju_df.columns else None
        name_col = 'taxon_name' if 'taxon_name' in kaiju_df.columns else None
        if pct_col and name_col:
            top = kaiju_df.nlargest(top_n, pct_col)
            for _, row in top.iterrows():
                name = str(row[name_col]).strip()
                if name.lower() in ('unclassified', 'cannot be assigned to a (non-viral) taxon',
                                    'cannot be assigned', 'na', 'nan'):
                    continue
                data.append({
                    'Species': name[:40],
                    'Tool': 'Kaiju',
                    'Abundance': float(row[pct_col])
                })

    if not data:
        return None

    df = pd.DataFrame(data)
    fig = go.Figure()
    tool_colors = {'Bracken': BLUE, 'Kraken2': ORANGE, 'Kaiju': GREEN}

    for tool in ['Bracken', 'Kraken2', 'Kaiju']:
        sub = df[df['Tool'] == tool].sort_values('Abundance', ascending=True)
        if not sub.empty:
            fig.add_trace(go.Bar(
                name=tool, y=sub['Species'], x=sub['Abundance'],
                orientation='h', marker_color=tool_colors[tool], opacity=0.85,
                text=[f'{v:.1f}%' for v in sub['Abundance']],
                textposition='outside', textfont_size=8
            ))

    n_bars = min(len(df), top_n * 3)
    fig.update_layout(
        title=dict(text=f'Top species — {sample}', font_size=14, font_color=BLUE),
        barmode='group',
        height=max(300, n_bars * 22 + 120),
        margin=dict(l=220, r=80, t=50, b=40),
        xaxis_title='Relative abundance (%)',
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=10),
    )
    return fig


def chart_rarefaction(sylph_df, samples):
    """Curvas de rarefacción por muestra (estimadas desde Sylph GTDB)."""
    fig = go.Figure()
    bact = sylph_df[sylph_df['DB'] == 'Bacteria/Archaea']

    any_data = False
    for i, sample in enumerate(sorted(samples)):
        s = bact[bact['Sample'] == sample]
        if s.empty or len(s) < 2:
            continue

        grp = s.groupby('Species')['Tax_abund'].sum().reset_index()
        grp = grp[grp['Tax_abund'] > 0]
        if len(grp) < 2:
            continue

        steps, means, lows, highs, n_total = simulate_rarefaction(
            grp['Tax_abund'].values, seed=42 + i
        )
        color = PALETTE[i % len(PALETTE)]

        # Banda de confianza
        r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
        fig.add_trace(go.Scatter(
            x=np.concatenate([steps, steps[::-1]]),
            y=np.concatenate([highs, lows[::-1]]),
            fill='toself',
            fillcolor=f'rgba({r},{g},{b},0.1)',
            line=dict(width=0), showlegend=False, hoverinfo='skip'
        ))
        # Línea media
        fig.add_trace(go.Scatter(
            x=steps, y=means, mode='lines',
            name=f'{sample} ({n_total} spp.)',
            line=dict(color=color, width=2.5),
        ))
        any_data = True

    if not any_data:
        return None

    fig.update_layout(
        title=dict(text='Rarefaction Curves (simulated from Sylph GTDB)',
                   font_size=14, font_color=BLUE),
        xaxis_title='Simulated reads',
        yaxis_title='Species detected',
        xaxis_type='log',
        height=500,
        margin=dict(l=60, r=30, t=60, b=60),
        legend=dict(font_size=9),
        font=dict(family='Arial, sans-serif'),
    )
    return fig


def chart_mag_viability(viability_df, min_coverage):
    """Gráfico de cobertura efectiva por organismo, con umbral de MAG."""
    if viability_df.empty:
        return None

    show = viability_df[viability_df['Eff_cov_sylph'] >= 5].copy()
    if show.empty:
        show = viability_df.nlargest(10, 'Eff_cov_sylph')
    if show.empty:
        return None

    fig = go.Figure()
    for i, sample in enumerate(sorted(show['Sample'].unique())):
        sub = show[show['Sample'] == sample].nlargest(8, 'Eff_cov_sylph')
        color = PALETTE[i % len(PALETTE)]
        fig.add_trace(go.Bar(
            name=sample,
            y=[f"{row['Species'][:30]} ({sample})" for _, row in sub.iterrows()],
            x=sub['Eff_cov_sylph'],
            orientation='h', marker_color=color, opacity=0.85,
            text=[f'{v:.1f}x' for v in sub['Eff_cov_sylph']],
            textposition='outside', textfont_size=9
        ))

    fig.add_vline(
        x=min_coverage, line_dash='dash', line_color=RED, line_width=2,
        annotation_text=f'Threshold {min_coverage}x',
        annotation_position='top right', annotation_font_color=RED
    )

    n_rows = len(show.drop_duplicates(['Sample', 'Species']))
    fig.update_layout(
        title=dict(
            text=f'Sylph Effective Coverage — MAG threshold: {min_coverage}x',
            font_size=14, font_color=BLUE),
        xaxis_title='Effective coverage (x)',
        height=max(400, n_rows * 28 + 150),
        margin=dict(l=280, r=80, t=60, b=60),
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=10),
    )
    return fig


def chart_fastqscreen(fqs_df):
    """Heatmap: % mapped per genome across samples."""
    if fqs_df.empty:
        return None

    pivot = fqs_df.pivot_table(index='Sample', columns='Genome',
                               values='Pct_mapped', fill_value=0)
    # Solo genomas con >0.05% en al menos una muestra
    active = pivot.columns[pivot.max() > 0.05].tolist()
    if not active:
        return None
    pivot = pivot[active]
    pivot = pivot.reindex(sorted(pivot.index))
    pivot = pivot[sorted(active, key=lambda g: -pivot[g].max())]

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale=[
            [0, '#FFFFFF'], [0.001, '#F0F4F7'], [0.01, '#B3D4E8'],
            [0.05, '#5B9ABF'], [0.2, '#2C5F8A'], [1.0, '#1A3F5C']
        ],
        text=[[f'{v:.2f}%' if v > 0 else '' for v in row] for row in pivot.values],
        texttemplate='%{text}',
        textfont_size=7,
        hovertemplate='%{y}<br>%{x}: %{z:.3f}%<extra></extra>',
        colorbar=dict(title=dict(text='% Mapped', side='right')),
        zmin=0, zmax=max(pivot.values.max(), 1.0)
    ))

    fig.update_layout(
        title=dict(text='FastQ Screen — Reads Mapped per Reference Genome',
                   font_size=14, font_color=BLUE),
        xaxis=dict(tickangle=-45, tickfont_size=7, side='bottom'),
        yaxis=dict(tickfont_size=9),
        height=max(350, len(pivot) * 28 + 200),
        margin=dict(l=160, r=100, t=60, b=180),
        font=dict(family='Arial, sans-serif'),
    )
    return fig


def chart_kma_amr(kma_df, min_identity=80, min_coverage=60):
    """Horizontal bar chart: resistance genes detected per sample."""
    if kma_df.empty:
        return None

    filt = kma_df[(kma_df['Template_Identity'] >= min_identity) &
                  (kma_df['Template_Coverage'] >= min_coverage)].copy()
    if filt.empty:
        return None

    # Top 20 genes por profundidad global
    top_genes = filt.groupby('Gene')['Depth'].max().nlargest(20).index.tolist()
    show = filt[filt['Gene'].isin(top_genes)].copy()

    fig = go.Figure()
    for i, sample in enumerate(sorted(show['Sample'].unique())):
        sub = show[show['Sample'] == sample].sort_values('Depth', ascending=True)
        color = PALETTE[i % len(PALETTE)]
        fig.add_trace(go.Bar(
            name=sample,
            y=sub['Gene'].str[:50], x=sub['Depth'],
            orientation='h', marker_color=color, opacity=0.85,
            text=[f'{v:.1f}x' for v in sub['Depth']],
            textposition='outside', textfont_size=8
        ))

    n_bars = len(show.drop_duplicates(['Sample', 'Gene']))
    fig.update_layout(
        title=dict(
            text=f'Resistance Genes (KMA ResFinder, id ≥{min_identity}%, cov ≥{min_coverage}%)',
            font_size=14, font_color=BLUE),
        xaxis_title='Depth (x)',
        barmode='group',
        height=max(400, n_bars * 22 + 150),
        margin=dict(l=350, r=80, t=60, b=60),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=10),
    )
    return fig


def classify_amr_class(gene_name):
    """Classify AMR gene into antibiotic class from gene name prefix."""
    g = str(gene_name).lower().split('_')[0]
    if 'bla' in g: return 'Beta-lactam'
    if any(x in g for x in ['aac', 'aad', 'ant', 'aph', 'str', 'rmt']): return 'Aminoglycoside'
    if any(x in g for x in ['erm', 'mef', 'msr', 'mph']): return 'Macrolide'
    if 'tet' in g: return 'Tetracycline'
    if any(x in g for x in ['qnr', 'oqx', 'qep']): return 'Quinolone'
    if 'cat' in g: return 'Phenicol'
    if 'sul' in g: return 'Sulfonamide'
    if 'dfr' in g: return 'Trimethoprim'
    if 'van' in g: return 'Glycopeptide'
    if 'mcr' in g: return 'Colistin'
    if 'fos' in g: return 'Fosfomycin'
    if any(x in g for x in ['lnu', 'lsa']): return 'Lincosamide'
    if any(x in g for x in ['cfr', 'optr', 'poxt']): return 'Oxazolidinone'
    return 'Other'


# Paleta fija para clases AMR (consistente entre runs)
AMR_CLASS_COLORS = {
    'Beta-lactam': '#E74C3C', 'Aminoglycoside': '#3498DB', 'Macrolide': '#2ECC71',
    'Tetracycline': '#F39C12', 'Quinolone': '#9B59B6', 'Phenicol': '#1ABC9C',
    'Sulfonamide': '#E67E22', 'Trimethoprim': '#34495E', 'Glycopeptide': '#C0392B',
    'Colistin': '#D35400', 'Fosfomycin': '#7F8C8D', 'Lincosamide': '#16A085',
    'Oxazolidinone': '#8E44AD', 'Other': '#95A5A6'
}


def chart_amr_classes(kma_df, min_identity=80, min_coverage=60):
    """Stacked bar chart: AMR determinants per antibiotic class per sample."""
    if kma_df.empty:
        return None

    filt = kma_df[(kma_df['Template_Identity'] >= min_identity) &
                  (kma_df['Template_Coverage'] >= min_coverage)].copy()
    if filt.empty:
        return None

    filt['Class'] = filt['Gene'].apply(classify_amr_class)
    pivot = filt.groupby(['Sample', 'Class']).size().unstack(fill_value=0)

    # Ordenar clases por total descendente
    class_order = pivot.sum().sort_values(ascending=False).index.tolist()
    samples = sorted(pivot.index.tolist())

    fig = go.Figure()
    for cls in class_order:
        if cls not in pivot.columns:
            continue
        color = AMR_CLASS_COLORS.get(cls, '#95A5A6')
        fig.add_trace(go.Bar(
            name=cls, x=samples, y=pivot.loc[samples, cls].tolist(),
            marker_color=color,
            hovertemplate=f'<b>{cls}</b><br>%{{x}}<br>%{{y}} genes<extra></extra>',
        ))

    fig.update_layout(
        title=dict(text='AMR Determinants by Antibiotic Class',
                   font_size=14, font_color=BLUE),
        xaxis_title='Sample', yaxis_title='No. resistance genes',
        barmode='stack',
        height=450,
        margin=dict(l=60, r=40, t=60, b=100),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5,
                    font_size=9),
        font=dict(family='Arial, sans-serif', size=11),
    )
    return fig


def chart_whittaker(sylph_df, samples):
    """Curvas rank-abundance (Whittaker) desde datos Sylph."""
    if sylph_df.empty:
        return None
    bact = sylph_df[sylph_df['DB'] == 'Bacteria/Archaea']
    fig = go.Figure()
    for i, sample in enumerate(samples):
        s_df = bact[bact['Sample'] == sample].copy()
        if s_df.empty:
            continue
        grp = s_df.groupby('Species')['Tax_abund'].sum().sort_values(ascending=False).reset_index()
        grp = grp[grp['Tax_abund'] > 0]
        if grp.empty:
            continue
        ranks = list(range(1, len(grp) + 1))
        fig.add_trace(go.Scatter(
            x=ranks, y=grp['Tax_abund'].tolist(),
            mode='lines+markers', name=sample,
            text=grp['Species'].tolist(),
            hovertemplate='<b>%{text}</b><br>Rank: %{x}<br>Abundance: %{y:.3f}%<extra>%{fullData.name}</extra>',
            marker=dict(size=4),
            line=dict(color=PALETTE[i % len(PALETTE)])
        ))
    fig.update_layout(
        title=dict(text='Rank-Abundance Curves (Whittaker) — Sylph GTDB',
                   font_size=14, font_color=BLUE),
        xaxis_title='Species rank',
        yaxis_title='Relative abundance (%)',
        yaxis_type='log',
        height=480,
        margin=dict(l=80, r=40, t=60, b=60),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=11),
    )
    return fig


def chart_depth_candidates(viability_df, min_coverage):
    """Dot-plot de cobertura estimada por especie × muestra con umbral de ensamblaje."""
    if viability_df.empty:
        return None
    # Filtrar a >= 1x para reducir ruido
    df = viability_df[viability_df['Est_cov'] >= 1].copy()
    if df.empty:
        return None
    # Tomar top N por muestra para no saturar
    parts = []
    for sample in df['Sample'].unique():
        parts.append(df[df['Sample'] == sample].nlargest(20, 'Est_cov'))
    df = pd.concat(parts).drop_duplicates()
    samples = sorted(df['Sample'].unique())
    fig = go.Figure()
    for i, sample in enumerate(samples):
        s_df = df[df['Sample'] == sample]
        fig.add_trace(go.Scatter(
            x=s_df['Est_cov'], y=s_df['Species'],
            mode='markers', name=sample,
            marker=dict(size=10, color=PALETTE[i % len(PALETTE)],
                        line=dict(width=1, color='white')),
            hovertemplate=(
                '<b>%{y}</b><br>%{fullData.name}<br>'
                'Est. coverage: %{x:.1f}x<br>'
                '<extra></extra>'
            ),
        ))
    # Linea umbral 30x
    fig.add_vline(x=min_coverage, line_dash='dash', line_color=RED, line_width=2,
                  annotation_text=f'{min_coverage:.0f}x', annotation_position='top right',
                  annotation_font_color=RED)
    species_list = sorted(df['Species'].unique())
    fig.update_layout(
        title=dict(text=f'MAG Assembly Candidates (threshold {min_coverage:.0f}x)',
                   font_size=14, font_color=BLUE),
        xaxis_title='Estimated coverage (x)', xaxis_type='log',
        height=max(400, 28 * len(species_list) + 120),
        margin=dict(l=250, r=60, t=60, b=60),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
        font=dict(family='Arial, sans-serif', size=11),
    )
    return fig


def chart_phenotypic(pheno_df):
    """Heatmap de breadth ≥1x por organismo × muestra con minimap2."""
    if pheno_df.empty:
        return None

    samples = sorted(pheno_df['sample'].unique())
    organisms = sorted(pheno_df['species'].unique())

    z, hover = [], []
    for org in organisms:
        row_z, row_h = [], []
        for sample in samples:
            mask = (pheno_df['sample'] == sample) & (pheno_df['species'] == org)
            if mask.any():
                r = pheno_df[mask].iloc[0]
                row_z.append(r['breadth_1x'])
                row_h.append(
                    f"<b>{org}</b><br>{sample}<br>"
                    f"Breadth ≥1x: {r['breadth_1x']:.1f}%<br>"
                    f"Breadth ≥10x: {r['breadth_10x']:.1f}%<br>"
                    f"Breadth ≥30x: {r['breadth_30x']:.1f}%<br>"
                    f"Mean depth: {r['mean_depth']:.1f}x<br>"
                    f"Mapped reads: {int(r['mapped']):,}"
                )
            else:
                row_z.append(0)
                row_h.append('No data')
        z.append(row_z)
        hover.append(row_h)

    # Text annotations: breadth_1x inside cells
    text = []
    for row in z:
        text.append([f'{v:.0f}%' if v > 0 else '' for v in row])

    fig = go.Figure(data=go.Heatmap(
        z=z, x=samples, y=organisms,
        text=text, texttemplate='%{text}', textfont=dict(size=11),
        customdata=hover,
        hovertemplate='%{customdata}<extra></extra>',
        colorscale=[
            [0.0, '#FFFFFF'], [0.05, '#FFEBEE'], [0.10, '#FFCDD2'],
            [0.25, '#FFE082'], [0.50, '#66BB6A'], [0.75, '#2E7D32'],
            [1.0, '#1B5E20'],
        ],
        zmin=0, zmax=100,
        colorbar=dict(title=dict(text='Breadth ≥1x (%)', side='right'))
    ))

    fig.update_layout(
        title=dict(text='Phenotypic Verification — minimap2 Coverage vs Type-Strains',
                   font_size=14, font_color=BLUE),
        xaxis_title='Sample', yaxis_title='Organism',
        height=max(350, 55 * len(organisms) + 120),
        margin=dict(l=200, r=100, t=60, b=60),
        font=dict(family='Arial, sans-serif', size=11),
    )
    return fig


def build_phenotypic_table(pheno_df):
    """Tabla HTML con métricas de cobertura fenotípica por muestra × organismo."""
    if pheno_df.empty:
        return '<p class="no-data">No target organisms provided (--phenotypic_targets).</p>'

    def breadth_color(val):
        """Color de fondo según breadth of coverage."""
        if val >= 70: return 'background:#2E7D32; color:white'
        if val >= 30: return 'background:#66BB6A; color:white'
        if val >= 10: return 'background:#FFE082'
        if val >= 1:  return 'background:#FFCDD2'
        return ''

    rows_html = []
    for _, r in pheno_df.sort_values(['sample', 'species']).iterrows():
        b1  = r['breadth_1x']
        b10 = r['breadth_10x']
        b30 = r['breadth_30x']
        rows_html.append(
            f'<tr>'
            f'<td>{r["sample"]}</td>'
            f'<td><i>{r["species"]}</i></td>'
            f'<td style="text-align:right">{int(r["mapped"]):,}</td>'
            f'<td style="text-align:right">{r["mean_depth"]:.1f}x</td>'
            f'<td style="text-align:right; {breadth_color(b1)}">{b1:.1f}%</td>'
            f'<td style="text-align:right; {breadth_color(b10)}">{b10:.1f}%</td>'
            f'<td style="text-align:right; {breadth_color(b30)}">{b30:.1f}%</td>'
            f'</tr>'
        )

    return f"""
    <table class="data-table">
        <thead><tr>
            <th>Sample</th><th>Organism</th><th>Mapped Reads</th>
            <th>Mean Depth</th><th>Breadth &ge;1x</th>
            <th>Breadth &ge;10x</th><th>Breadth &ge;30x</th>
        </tr></thead>
        <tbody>{''.join(rows_html)}</tbody>
    </table>
    <p style="font-size:0.85rem; color:#666; margin-top:0.5rem;">
        <b>Breadth</b> = percentage of the reference genome (type-strain) covered at the indicated
        depth threshold. Dark green &ge;70%, green &ge;30%, yellow &ge;10%,
        red &lt;10%. Mean depth is the average across all genome positions.
    </p>"""


# ════════════════════════════════════════════════════════════════════
# TABLAS HTML
# ════════════════════════════════════════════════════════════════════

def fig_to_html(fig):
    """Convierte figura Plotly a div HTML (sin JS incluido)."""
    if fig is None:
        return '<p class="no-data">No data available for this section.</p>'
    return fig.to_html(full_html=False, include_plotlyjs=False)


def build_survival_table(survival_df):
    """Tabla HTML de métricas de supervivencia."""
    rows_html = ''
    for _, r in survival_df.iterrows():
        pct_cls = 'good' if r['pct_reads'] >= 70 else ('warn' if r['pct_reads'] >= 50 else 'bad')
        gb_cls = 'highlight' if r['clean_gb'] >= 1.0 else ''
        rows_html += (
            f'<tr>'
            f'<td class="sample">{r["Sample"]}</td>'
            f'<td>{r["raw_reads"]:,.0f}</td>'
            f'<td>{r["raw_gb"]:.2f}</td>'
            f'<td>{r["clean_reads"]:,.0f}</td>'
            f'<td class="{gb_cls}">{r["clean_gb"]:.2f}</td>'
            f'<td class="{pct_cls}">{r["pct_reads"]:.1f}%</td>'
            f'<td class="{pct_cls}">{r["pct_gb"]:.1f}%</td>'
            f'</tr>\n'
        )

    # Fila de totales
    tot = survival_df[['raw_reads', 'raw_gb', 'clean_reads', 'clean_gb']].sum()
    t_pct_r = (tot['clean_reads'] / tot['raw_reads'] * 100) if tot['raw_reads'] > 0 else 0
    t_pct_g = (tot['clean_gb'] / tot['raw_gb'] * 100) if tot['raw_gb'] > 0 else 0
    rows_html += (
        f'<tr class="total-row">'
        f'<td><strong>TOTAL</strong></td>'
        f'<td><strong>{tot["raw_reads"]:,.0f}</strong></td>'
        f'<td><strong>{tot["raw_gb"]:.2f}</strong></td>'
        f'<td><strong>{tot["clean_reads"]:,.0f}</strong></td>'
        f'<td><strong>{tot["clean_gb"]:.2f}</strong></td>'
        f'<td><strong>{t_pct_r:.1f}%</strong></td>'
        f'<td><strong>{t_pct_g:.1f}%</strong></td>'
        f'</tr>\n'
    )

    return (
        '<table class="data-table">'
        '<thead><tr>'
        '<th>Sample</th><th>Raw Reads</th><th>Raw Gb</th>'
        '<th>Clean Reads</th><th>Clean Gb</th>'
        '<th>% Reads Ret.</th><th>% Gb Ret.</th>'
        '</tr></thead>'
        f'<tbody>{rows_html}</tbody>'
        '</table>'
    )


def build_viability_table(viability_df, min_coverage):
    """Tabla HTML de viabilidad para MAGs."""
    if viability_df.empty:
        return '<p class="no-data">No viability data available.</p>'

    rows_html = ''
    for sample in sorted(viability_df['Sample'].unique()):
        sub = viability_df[viability_df['Sample'] == sample].nlargest(10, 'Eff_cov_sylph')
        for _, r in sub.iterrows():
            if r['Eff_cov_sylph'] >= min_coverage:
                badge = '<span class="badge badge-green">VIABLE</span>'
            elif r['Eff_cov_sylph'] >= 10:
                badge = '<span class="badge badge-orange">POSSIBLE</span>'
            elif r['Eff_cov_sylph'] >= 5:
                badge = '<span class="badge badge-grey">MARGINAL</span>'
            else:
                continue
            rows_html += (
                f'<tr>'
                f'<td class="sample">{r["Sample"]}</td>'
                f'<td><em>{r["Species"]}</em></td>'
                f'<td>{r["Abund_pct"]:.1f}%</td>'
                f'<td><strong>{r["Eff_cov_sylph"]:.1f}x</strong></td>'
                f'<td>{r["Est_cov"]:.1f}x</td>'
                f'<td>{badge}</td>'
                f'</tr>\n'
            )

    if not rows_html:
        return '<p class="no-data">No organism exceeds 5x coverage.</p>'

    return (
        '<table class="data-table">'
        '<thead><tr>'
        '<th>Sample</th><th>Species</th><th>Abundance</th>'
        f'<th>Sylph Cov. (x)</th><th>Est. Cov. (x)</th><th>MAG {min_coverage}x</th>'
        '</tr></thead>'
        f'<tbody>{rows_html}</tbody>'
        '</table>'
        f'<p class="note">'
        f'Sylph Cov. = effective coverage reported directly by Sylph. '
        f'Est. Cov. = Clean Gb &times; abundance / average genome size. '
        f'VIABLE &ge;{min_coverage}x | POSSIBLE &ge;10x | MARGINAL &ge;5x'
        f'</p>'
    )


def build_conclusion(survival_df, viability_df, min_coverage):
    """Genera texto dinámico de conclusión logística."""
    n_samples = len(survival_df)
    n_above_1gb = int((survival_df['clean_gb'] >= 1.0).sum())
    total_clean_gb = survival_df['clean_gb'].sum()
    mean_retention = survival_df['pct_reads'].mean()

    viable = viability_df[viability_df['Eff_cov_sylph'] >= min_coverage] if not viability_df.empty else pd.DataFrame()
    n_viable = len(viable)
    n_viable_samples = viable['Sample'].nunique() if not viable.empty else 0

    parts = []

    parts.append(
        f'Of the <strong>{n_samples} samples</strong> processed, '
        f'<strong>{n_above_1gb}</strong> exceed 1 Gb of clean data, '
        f'making them suitable for MAG recovery.'
    )
    parts.append(
        f'Total clean data volume is <strong>{total_clean_gb:.2f} Gb</strong> '
        f'with a mean retention of <strong>{mean_retention:.1f}%</strong> after '
        f'Porechop ABI + Chopper.'
    )

    if n_viable > 0:
        parts.append(
            f'<strong>{n_viable} organisms</strong> were identified in '
            f'<strong>{n_viable_samples} samples</strong> with coverage '
            f'&ge;{min_coverage}x according to Sylph, viable for targeted assembly.'
        )
    else:
        parts.append(
            f'No organism reaches the {min_coverage}x threshold for direct assembly. '
            f'Consider increasing sequencing depth (fewer samples per flow cell).'
        )

    if n_above_1gb >= n_samples * 0.7:
        parts.append(
            '<span class="rec-good">Recommendation: the current sample density per flow cell '
            'is adequate for taxonomic profiling and recovery of dominant MAGs.</span>'
        )
    elif n_above_1gb >= n_samples * 0.4:
        parts.append(
            '<span class="rec-warn">Recommendation: consider reducing to 5-7 samples per flow cell '
            'to improve individual depth and MAG viability.</span>'
        )
    else:
        parts.append(
            '<span class="rec-bad">Recommendation: significantly reduce the number of samples '
            'per flow cell (3-5 maximum) to achieve sufficient depth for MAGs.</span>'
        )

    return '\n'.join(f'<p>{p}</p>' for p in parts)


def build_fastqscreen_table(fqs_df):
    """Tabla resumen de FastQ Screen: genomas con mayor % de mapping."""
    if fqs_df.empty:
        return '<p class="no-data">No FastQ Screen data available.</p>'

    rows_html = ''
    for sample in sorted(fqs_df['Sample'].unique()):
        sub = fqs_df[fqs_df['Sample'] == sample].nlargest(5, 'Pct_mapped')
        for _, r in sub.iterrows():
            if r['Pct_mapped'] < 0.01:
                continue
            cls = 'warn' if r['Pct_mapped'] >= 1.0 else ''
            rows_html += (
                f'<tr>'
                f'<td class="sample">{r["Sample"]}</td>'
                f'<td>{r["Genome"]}</td>'
                f'<td class="{cls}"><strong>{r["Pct_mapped"]:.2f}%</strong></td>'
                f'<td>{r["Pct_unique"]:.2f}%</td>'
                f'<td>{r["Pct_multi"]:.2f}%</td>'
                f'</tr>\n'
            )

    if not rows_html:
        return '<p class="no-data">No reference genome detected.</p>'

    return (
        '<table class="data-table">'
        '<thead><tr>'
        '<th>Sample</th><th>Genome</th><th>% Total Mapped</th>'
        '<th>% Unique Hit</th><th>% Multi-hit</th>'
        '</tr></thead>'
        f'<tbody>{rows_html}</tbody>'
        '</table>'
        '<p class="note">Top 5 genomes per sample. '
        'Values &ge;1% are highlighted. Screening performed on filtered reads.</p>'
    )


def build_kma_table(kma_df, min_identity=80, min_coverage=60):
    """Tabla de genes de resistencia detectados por KMA."""
    if kma_df.empty:
        return '<p class="no-data">No AMR resistance data available.</p>'

    filt = kma_df[(kma_df['Template_Identity'] >= min_identity) &
                  (kma_df['Template_Coverage'] >= min_coverage)].copy()

    if filt.empty:
        return (f'<p class="no-data">No resistance gene detected with '
                f'identity &ge;{min_identity}% and coverage &ge;{min_coverage}%.</p>')

    filt = filt.sort_values(['Sample', 'Depth'], ascending=[True, False])

    rows_html = ''
    for _, r in filt.iterrows():
        depth_cls = 'good' if r['Depth'] >= 10 else ('warn' if r['Depth'] >= 1 else '')
        rows_html += (
            f'<tr>'
            f'<td class="sample">{r["Sample"]}</td>'
            f'<td><em>{r["Gene"][:60]}</em></td>'
            f'<td>{r["Template_Identity"]:.1f}%</td>'
            f'<td>{r["Template_Coverage"]:.1f}%</td>'
            f'<td class="{depth_cls}"><strong>{r["Depth"]:.1f}x</strong></td>'
            f'</tr>\n'
        )

    return (
        '<table class="data-table">'
        '<thead><tr>'
        '<th>Sample</th><th>Resistance Gene</th><th>Identity (%)</th>'
        '<th>Coverage (%)</th><th>Depth (x)</th>'
        '</tr></thead>'
        f'<tbody>{rows_html}</tbody>'
        '</table>'
        f'<p class="note">'
        f'Filter: identity &ge;{min_identity}%, coverage &ge;{min_coverage}%. '
        f'Depth &ge;10x indicates high confidence. '
        f'DB: ResFinder (DTU). Direct screening from reads (KMA -ont).'
        f'</p>'
    )


# ════════════════════════════════════════════════════════════════════
# ENSAMBLAJE HTML
# ════════════════════════════════════════════════════════════════════

# CSS como constante (no f-string) para evitar conflictos con {}
HTML_CSS = """<style>
:root {
    --blue: #2C5F8A; --blue-light: #5B9ABF; --orange: #E05C2A;
    --green: #28A745; --red: #E74C3C; --grey: #6C757D;
    --bg: #F8F9FA; --card-bg: #FFFFFF;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
    font-family: 'Segoe UI', Tahoma, Arial, sans-serif;
    background: var(--bg); color: #333; line-height: 1.6;
}
.header {
    background: linear-gradient(135deg, var(--blue) 0%, #1a3f5c 100%);
    color: white; padding: 2rem 3rem;
    border-bottom: 4px solid var(--orange);
}
.header h1 { font-size: 1.8rem; margin-bottom: 0.3rem; }
.header .subtitle { opacity: 0.8; font-size: 0.95rem; }
.container { max-width: 1200px; margin: 0 auto; padding: 1.5rem; }
.section {
    background: var(--card-bg); border-radius: 8px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    margin-bottom: 1.5rem; overflow: hidden;
}
.section-header {
    background: var(--blue); color: white;
    padding: 0.8rem 1.5rem; font-size: 1.1rem; font-weight: 600;
}
.section-body { padding: 1.5rem; }
.section-body > p { margin-bottom: 0.8rem; color: #555; font-size: 0.92rem; }
.data-table {
    width: 100%; border-collapse: collapse;
    font-size: 0.85rem; margin: 1rem 0;
}
.data-table th {
    background: var(--blue); color: white;
    padding: 0.6rem 0.8rem; text-align: center;
    font-weight: 600; font-size: 0.8rem;
}
.data-table td {
    padding: 0.5rem 0.8rem; text-align: center;
    border-bottom: 1px solid #eee;
}
.data-table tr:nth-child(even) { background: #F0F4F7; }
.data-table tr:hover { background: #E8EFF5; }
.data-table .sample {
    text-align: left; font-weight: 600;
    font-family: 'Courier New', monospace; font-size: 0.82rem;
}
.data-table .good { color: var(--green); font-weight: 600; }
.data-table .warn { color: var(--orange); font-weight: 600; }
.data-table .bad { color: var(--red); font-weight: 600; }
.data-table .highlight { background: #C8E6C9 !important; font-weight: 700; }
.data-table .total-row { background: #E3EBF3 !important; font-size: 0.9rem; }
.badge {
    display: inline-block; padding: 0.2rem 0.6rem;
    border-radius: 12px; font-size: 0.75rem;
    font-weight: 700; color: white;
}
.badge-green { background: var(--green); }
.badge-orange { background: var(--orange); }
.badge-grey { background: var(--grey); }
.note { font-size: 0.8rem; color: var(--grey); margin-top: 0.8rem; }
.no-data { color: #999; font-style: italic; padding: 1rem 0; }
.conclusion {
    background: #E8F5E9; border-left: 4px solid var(--green);
    padding: 1.2rem; border-radius: 4px; margin-top: 0.5rem;
}
.conclusion p { margin-bottom: 0.6rem; }
.rec-good { color: var(--green); font-weight: 700; }
.rec-warn { color: var(--orange); font-weight: 700; }
.rec-bad { color: var(--red); font-weight: 700; }
.tax-separator {
    border: none; border-top: 1px solid #eee;
    margin: 1.5rem 0;
}
.footer {
    text-align: center; padding: 1.5rem;
    color: var(--grey); font-size: 0.8rem;
}
.plotly-graph-div { margin: 0 auto; }
.sample-filter {
    background: #fff; border-radius: 8px; padding: 0.8rem 1.2rem;
    margin-bottom: 1.2rem; box-shadow: 0 2px 6px rgba(0,0,0,0.08);
    display: flex; flex-wrap: wrap; align-items: center; gap: 0.6rem;
}
.sample-filter strong { color: var(--blue); font-size: 0.9rem; margin-right: 0.5rem; }
.sample-filter label {
    font-size: 0.82rem; cursor: pointer; padding: 0.2rem 0.5rem;
    border-radius: 4px; transition: background 0.15s;
}
.sample-filter label:hover { background: #E8EFF5; }
.data-table th { cursor: pointer; user-select: none; }
.data-table th:after { content: ' \\21C5'; opacity: 0.3; font-size: 0.7rem; }
.data-table th.sort-asc:after { content: ' \\25B2'; opacity: 0.8; }
.data-table th.sort-desc:after { content: ' \\25BC'; opacity: 0.8; }
@media print {
    .section { break-inside: avoid; }
    .sample-filter { display: none; }
    body { background: white; }
}
</style>"""


def build_html(sections, run_name, plotly_js):
    """Ensambla el reporte HTML completo y autocontenido."""
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')
    min_q = sections['min_quality']
    min_l = sections['min_length']
    gs = sections['genome_size']
    mc = sections['min_coverage']

    # Construir body con f-string (CSS y plotly_js son variables, sus {} internos son seguros)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Metagenomic Report — {run_name}</title>
    <script type="text/javascript">{plotly_js}</script>
    {HTML_CSS}
</head>
<body>
    <div class="header">
        <h1>EpiTax — Nanopore Metagenomic Report</h1>
        <div class="subtitle">Run: {run_name} &nbsp;|&nbsp; Generated: {date_str} &nbsp;|&nbsp; FISABIO — EPIMOL</div>
    </div>
    <div class="container">

        <!-- Sample filter -->
        <div id="sample-filter" class="sample-filter">
            <strong>Filter samples:</strong>
        </div>

        <!-- 1. Survival -->
        <div class="section">
            <div class="section-header">1. Survival Metrics (Reads and Bases)</div>
            <div class="section-body">
                <p>Comparison before (raw) and after (clean) filtering with Porechop ABI
                (adapter trimming) + Chopper (Q&ge;{min_q}, len&ge;{min_l}bp).
                Green cells indicate &ge;1 Gb of clean data.</p>
                {sections['survival_table']}
                {sections['survival_chart']}
            </div>
        </div>

        <!-- 2. FastQ Screen -->
        <div class="section">
            <div class="section-header">2. FastQ Screen — Contamination Screening</div>
            <div class="section-body">
                <p>Percentage of filtered reads mapping against reference genomes
                (host, controls, pathogens, fecal indicators, viruses, parasites).
                Values &ge;1% are potentially significant.</p>
                {sections['fastqscreen_table']}
                {sections['fastqscreen_chart']}
            </div>
        </div>

        <!-- 3. Integrated Taxonomy -->
        <div class="section">
            <div class="section-header">3. Integrated Taxonomic Summary</div>
            <div class="section-body">
                <p>Top species detected by <strong>Bracken</strong> (refined abundances
                from Kraken2), <strong>Kraken2</strong> (k-mers) and <strong>Kaiju</strong>
                (protein-level). Each chart shows an individual sample.</p>
                {sections['taxonomy_charts']}
            </div>
        </div>

        <!-- 4. Diversity — Whittaker -->
        <div class="section">
            <div class="section-header">4. Diversity Profile — Rank-Abundance Curves</div>
            <div class="section-body">
                <p><strong>Whittaker</strong> rank-abundance curves based on Sylph (GTDB).
                The X-axis shows species ranked from most to least abundant;
                the Y-axis shows their relative abundance (log scale).
                Longer curves = more species; flatter curves = higher evenness.</p>
                {sections['whittaker_chart']}
            </div>
        </div>

        <!-- 5. MAG Viability -->
        <div class="section">
            <div class="section-header">5. MAG Viability Estimation</div>
            <div class="section-body">
                <p>Organisms with sufficient coverage for genome assembly (MAGs).
                <strong>Sylph coverage</strong> is the direct estimate from the ANI profiler.
                <strong>Estimated coverage</strong> is calculated as:
                Clean Gb &times; abundance / genome size ({gs} Mb).
                Minimum threshold: <strong>{mc}x</strong>.
                The lower chart shows species with &ge;1x and the assembly
                threshold ({mc}x) as a red line.</p>
                {sections['viability_table']}
                {sections['viability_chart']}
                {sections.get('depth_chart', '')}
            </div>
        </div>

        <!-- 6. AMR Screening -->
        <div class="section">
            <div class="section-header">6. Antimicrobial Resistance (AMR) Screening</div>
            <div class="section-body">
                <p>Direct detection of antibiotic resistance genes from filtered reads
                (assembly-free). <strong>KMA</strong> (-ont mode for Nanopore)
                against the <strong>ResFinder</strong> database (DTU).</p>
                {sections.get('amr_class_chart', '')}
                {sections['kma_table']}
                {sections['kma_chart']}
            </div>
        </div>

        <!-- 7. Phenotypic Verification -->
        {sections.get('phenotypic_section', '')}

        <!-- 8. Rarefaction -->
        <div class="section">
            <div class="section-header">8. Rarefaction Curves</div>
            <div class="section-body">
                <p>Diversity saturation estimation based on abundances
                detected by Sylph (GTDB). If the curve flattens, the sequencing depth
                captures most of the diversity.
                Bands = P10-P90 interval from 10 simulations.</p>
                {sections['rarefaction_chart']}
            </div>
        </div>

        <!-- 9. Conclusion -->
        <div class="section">
            <div class="section-header">9. Conclusion and Recommendations</div>
            <div class="section-body">
                <div class="conclusion">
                    {sections['conclusion']}
                </div>
            </div>
        </div>

    </div>
    <div class="footer">
        EpiTax Pipeline v2 &nbsp;|&nbsp; FISABIO — EPIMOL &nbsp;|&nbsp; {date_str}
    </div>

    <script>
    // ── Tablas ordenables ──
    document.querySelectorAll('.data-table th').forEach(function(th) {{
        th.addEventListener('click', function() {{
            var table = th.closest('table');
            var idx = Array.from(th.parentNode.children).indexOf(th);
            var asc = !th.classList.contains('sort-asc');
            th.parentNode.querySelectorAll('th').forEach(function(h) {{
                h.classList.remove('sort-asc', 'sort-desc');
            }});
            th.classList.add(asc ? 'sort-asc' : 'sort-desc');
            var rows = Array.from(table.tBodies[0].rows);
            rows.sort(function(a, b) {{
                var x = a.cells[idx].textContent.trim();
                var y = b.cells[idx].textContent.trim();
                var xn = parseFloat(x.replace(/[,%x]/g, ''));
                var yn = parseFloat(y.replace(/[,%x]/g, ''));
                if (!isNaN(xn) && !isNaN(yn)) return asc ? xn - yn : yn - xn;
                return asc ? x.localeCompare(y) : y.localeCompare(x);
            }});
            rows.forEach(function(r) {{ table.tBodies[0].appendChild(r); }});
        }});
    }});

    // ── Filtro de muestras ──
    var allSamples = [];
    document.querySelectorAll('.data-table tbody tr').forEach(function(tr) {{
        var s = tr.cells[0] ? tr.cells[0].textContent.trim() : '';
        if (s && allSamples.indexOf(s) === -1) allSamples.push(s);
    }});
    allSamples.sort();
    var filterDiv = document.getElementById('sample-filter');
    if (filterDiv && allSamples.length > 1) {{
        allSamples.forEach(function(s) {{
            var lbl = document.createElement('label');
            var cb = document.createElement('input');
            cb.type = 'checkbox'; cb.checked = true; cb.value = s;
            cb.addEventListener('change', applySampleFilter);
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(' ' + s));
            filterDiv.appendChild(lbl);
        }});
    }} else if (filterDiv) {{ filterDiv.style.display = 'none'; }}

    function applySampleFilter() {{
        var checked = Array.from(document.querySelectorAll('#sample-filter input:checked'))
            .map(function(cb) {{ return cb.value; }});
        // Filtrar filas de tablas
        document.querySelectorAll('.data-table tbody tr').forEach(function(tr) {{
            if (tr.classList.contains('total-row')) return;
            var cell = tr.cells[0];
            if (!cell) return;
            var s = cell.textContent.trim();
            tr.style.display = checked.indexOf(s) >= 0 ? '' : 'none';
        }});
        // Filtrar gráficos por muestra
        document.querySelectorAll('[data-sample-chart]').forEach(function(el) {{
            var s = el.getAttribute('data-sample-chart');
            el.style.display = checked.indexOf(s) >= 0 ? '' : 'none';
        }});
    }}
    </script>
</body>
</html>"""


# ════════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Generate interactive HTML report for the Nanopore metagenomic pipeline.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--qc-raw', required=True, help='QC raw metrics CSV')
    parser.add_argument('--qc-filtered', required=True, help='QC filtered metrics CSV')
    parser.add_argument('--sylph', required=True, help='Sylph profile TSV (all samples)')
    parser.add_argument('--run-name', required=True, help='Run name')
    parser.add_argument('--output', required=True, help='Output HTML file')
    parser.add_argument('--genome-size', type=float, default=4.5,
                        help='Average bacterial genome size in Mb (default: 4.5)')
    parser.add_argument('--min-coverage', type=float, default=30,
                        help='Minimum coverage for MAGs (default: 30)')
    parser.add_argument('--min-quality', type=int, default=10,
                        help='Chopper quality threshold used (default: 10)')
    parser.add_argument('--min-length', type=int, default=1000,
                        help='Chopper length threshold used (default: 1000)')
    parser.add_argument('--phenotypic', default='',
                        help='Combined phenotypic verification TSV (minimap2, optional)')
    parser.add_argument('--results-dir', default='',
                        help='Results directory (reads from subdirectories; if empty, globs in cwd)')
    args = parser.parse_args()

    print(f'[INFO] Generating report for run: {args.run_name}')

    # ── Cargar datos ─────────────────────────────────────────────────
    qc_raw = parse_qc_csv(args.qc_raw)
    qc_filtered = parse_qc_csv(args.qc_filtered)
    sylph_df = parse_sylph_profile(args.sylph)
    samples = sorted(qc_filtered['Sample'].unique())
    print(f'[INFO] {len(samples)} samples detected')

    # Determinar directorios fuente (--results-dir o cwd)
    rd = args.results_dir
    kraken2_files = sorted(glob.glob(f'{rd}/07_tax_kraken2/*.kraken2.report') if rd
                           else glob.glob('*.kraken2.report'))
    bracken_files = sorted(glob.glob(f'{rd}/08_tax_bracken/*.bracken.S.txt') if rd
                           else glob.glob('*.bracken.S.txt'))
    kaiju_files = sorted(glob.glob(f'{rd}/09_tax_kaiju/*.kaiju.species.tsv') if rd
                         else glob.glob('*.kaiju.species.tsv'))
    fqs_files = sorted(glob.glob(f'{rd}/02_fastqscreen_raw/*_screen.txt') if rd
                       else glob.glob('*_screen.txt'))
    kma_res_files = sorted(glob.glob(f'{rd}/11_amr_kma/*.res') if rd
                           else glob.glob('*.res'))

    kraken2_data = {}
    for f in kraken2_files:
        sample = os.path.basename(f).replace('.kraken2.report', '')
        kraken2_data[sample] = parse_kraken2_report(f)

    bracken_data = {}
    for f in bracken_files:
        sample = os.path.basename(f).replace('.bracken.S.txt', '')
        bracken_data[sample] = parse_bracken_species(f)

    kaiju_data = {}
    for f in kaiju_files:
        sample = os.path.basename(f).replace('.kaiju.species.tsv', '')
        kaiju_data[sample] = parse_kaiju_species(f)

    print(f'[INFO] Kraken2: {len(kraken2_data)} samples | '
          f'Bracken: {len(bracken_data)} | Kaiju: {len(kaiju_data)} | Sylph: OK')

    fqs_df = parse_fastqscreen_results(fqs_files) if fqs_files else pd.DataFrame()
    print(f'[INFO] FastQ Screen: {len(fqs_files)} files')

    kma_df = parse_kma_results(kma_res_files) if kma_res_files else pd.DataFrame()
    n_kma_genes = len(kma_df[(kma_df['Template_Identity'] >= 80) &
                             (kma_df['Template_Coverage'] >= 60)]) if not kma_df.empty else 0
    print(f'[INFO] KMA: {len(kma_res_files)} files, {n_kma_genes} AMR genes (id>=80%, cov>=60%)')

    # Cargar verificación fenotípica (opcional)
    pheno_df = parse_phenotypic_results(args.phenotypic) if args.phenotypic else pd.DataFrame()
    if not pheno_df.empty:
        n_org = pheno_df['species'].nunique()
        print(f'[INFO] Phenotypic: {len(pheno_df)} entries, {n_org} organisms')
    else:
        print('[INFO] Phenotypic: not provided')

    # ── Análisis ─────────────────────────────────────────────────────
    survival = compute_survival(qc_raw, qc_filtered)
    viability = compute_mag_viability(
        sylph_df, survival, args.genome_size, args.min_coverage)

    # ── Gráficos ─────────────────────────────────────────────────────
    survival_chart = fig_to_html(chart_survival(survival))

    # FastQ Screen
    fastqscreen_chart = fig_to_html(chart_fastqscreen(fqs_df))

    # Taxonomía por muestra (envuelto en div filtrable)
    tax_parts = []
    for sample in samples:
        fig = chart_taxonomy_sample(
            sample,
            bracken_df=bracken_data.get(sample),
            kraken_df=kraken2_data.get(sample),
            kaiju_df=kaiju_data.get(sample),
        )
        if fig is not None:
            tax_parts.append(
                f'<div data-sample-chart="{sample}">{fig_to_html(fig)}</div>')
    taxonomy_html = ('\n<hr class="tax-separator">\n'.join(tax_parts)
                     if tax_parts
                     else '<p class="no-data">No taxonomic data available.</p>')

    # Whittaker rank-abundance
    whittaker_chart = fig_to_html(chart_whittaker(sylph_df, samples))

    rarefaction_html = fig_to_html(chart_rarefaction(sylph_df, samples))
    viability_chart = fig_to_html(chart_mag_viability(viability, args.min_coverage))

    # Candidatos MAG (dot-plot con umbral 30x)
    depth_chart = fig_to_html(chart_depth_candidates(viability, args.min_coverage))

    # KMA AMR
    kma_chart = fig_to_html(chart_kma_amr(kma_df))
    amr_class_chart = fig_to_html(chart_amr_classes(kma_df))

    # Verificación fenotípica (opcional)
    if not pheno_df.empty:
        phenotypic_chart = fig_to_html(chart_phenotypic(pheno_df))
        phenotypic_table = build_phenotypic_table(pheno_df)
        phenotypic_section = f"""
        <div class="section">
            <div class="section-header">7. Phenotypic Verification — minimap2 vs Type-Strains</div>
            <div class="section-body">
                <p>Direct mapping of filtered reads against reference genomes (type-strains, NCBI RefSeq)
                with <strong>minimap2</strong> (map-ont). The <b>breadth of coverage</b>
                (percentage of genome covered) is reported at three depth thresholds:
                &ge;1x (detection), &ge;10x (confidence), &ge;30x (assembly).
                Contigs from the same organism (chromosome + plasmids) are aggregated.</p>
                {phenotypic_table}
                {phenotypic_chart}
            </div>
        </div>"""
    else:
        phenotypic_section = ''

    # ── Plotly JS autocontenido ──────────────────────────────────────
    plotly_js = plotly.offline.get_plotlyjs()

    # ── Ensamblar HTML ───────────────────────────────────────────────
    sections = {
        'survival_table': build_survival_table(survival),
        'survival_chart': survival_chart,
        'fastqscreen_table': build_fastqscreen_table(fqs_df),
        'fastqscreen_chart': fastqscreen_chart,
        'taxonomy_charts': taxonomy_html,
        'whittaker_chart': whittaker_chart,
        'viability_table': build_viability_table(viability, args.min_coverage),
        'viability_chart': viability_chart,
        'depth_chart': depth_chart,
        'kma_table': build_kma_table(kma_df),
        'kma_chart': kma_chart,
        'amr_class_chart': amr_class_chart,
        'phenotypic_section': phenotypic_section,
        'rarefaction_chart': rarefaction_html,
        'conclusion': build_conclusion(survival, viability, args.min_coverage),
        'min_quality': args.min_quality,
        'min_length': args.min_length,
        'genome_size': args.genome_size,
        'min_coverage': args.min_coverage,
    }

    html = build_html(sections, args.run_name, plotly_js)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f'[OK] Report generated: {args.output}')


if __name__ == '__main__':
    main()
