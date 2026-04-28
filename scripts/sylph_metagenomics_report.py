#!/usr/bin/env python3
"""
sylph_metagenomics_report.py
Informe PDF de perfilado taxonómico Sylph — versión simplificada.

Centrado en GTDB (Bacteria/Archaea):
  - Composición taxonómica top especies
  - Curva Whittaker (rank-abundance)
  - Curva de rarefacción
  - Candidatos a ensamblaje (cobertura ≥ umbral)

Uso:
    python3 sylph_metagenomics_report.py <sylph_profile.tsv> [output.pdf] [--min-cov 5]
"""

import sys
import os
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

warnings.filterwarnings('ignore')

# ─── Paleta ──────────────────────────────────────────────────────────
PALETTE = [
    '#2C5F8A', '#5B9ABF', '#E05C2A', '#28A745', '#9B59B6',
    '#F39C12', '#1ABC9C', '#E74C3C', '#3498DB', '#95A5A6',
    '#D35400', '#2ECC71', '#8E44AD', '#16A085', '#C0392B',
]
GREY = '#CCCCCC'
BLUE = '#2C5F8A'
ORANGE = '#E05C2A'
GREEN = '#28A745'


# ─── Carga de datos ──────────────────────────────────────────────────

def extract_sample(path):
    base = os.path.basename(path)
    for ext in ['.filtered.fastq.gz', '.fastq.gz', '.fastq', '.fq.gz', '.fq']:
        base = base.replace(ext, '')
    return base


def extract_species(contig_name):
    """Extrae 'Genus species' del campo Contig_name de GTDB."""
    parts = contig_name.split()
    if len(parts) < 3:
        return contig_name[:40]
    skip = 1
    if len(parts) > skip and parts[skip].upper() == 'MAG:':
        skip += 1
    genus = parts[skip] if skip < len(parts) else ''
    species = parts[skip + 1] if skip + 1 < len(parts) else ''
    name = f'{genus} {species}'.strip()
    return name if name else contig_name[:40]


def classify_db(genome_file):
    gf = genome_file.lower()
    if 'imgvr' in gf:
        return 'Virus'
    elif 'fungi' in gf:
        return 'Fungi'
    else:
        return 'Bacteria/Archaea'


def load_data(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['Sample'] = df['Sample_file'].apply(extract_sample)
    df['DB'] = df['Genome_file'].apply(classify_db)
    df['Species'] = df['Contig_name'].apply(extract_species)
    df['Tax_abund'] = pd.to_numeric(df['Taxonomic_abundance'], errors='coerce').fillna(0)
    df['Eff_cov'] = pd.to_numeric(df['Eff_cov'], errors='coerce').fillna(0)
    df['Adj_ANI'] = pd.to_numeric(df['Adjusted_ANI'], errors='coerce').fillna(0)
    return df


def get_gtdb(df):
    """Filtra solo Bacteria/Archaea, agrupa por especie."""
    sub = df[df['DB'] == 'Bacteria/Archaea'].copy()
    return sub


# ─── Gráficos ────────────────────────────────────────────────────────

def fig_composition(sample_df, sample_name, top_n=15):
    """Barras horizontales de top especies (GTDB)."""
    gtdb = get_gtdb(sample_df)
    if gtdb.empty:
        return None

    grp = gtdb.groupby('Species').agg(
        abund=('Tax_abund', 'sum'),
        cov=('Eff_cov', 'mean'),
        ani=('Adj_ANI', 'mean')
    ).sort_values('abund', ascending=False).reset_index()

    top = grp.head(top_n).copy()
    rest = grp.iloc[top_n:]['abund'].sum()
    if rest > 0:
        top = pd.concat([top, pd.DataFrame([{
            'Species': f'Otros ({len(grp) - top_n} spp.)',
            'abund': rest, 'cov': 0, 'ani': 0
        }])], ignore_index=True)

    fig, ax = plt.subplots(figsize=(12, max(5, len(top) * 0.4)))
    colors = [PALETTE[i % len(PALETTE)] if i < top_n else GREY for i in range(len(top))]

    y = np.arange(len(top))
    bars = ax.barh(y, top['abund'], color=colors, edgecolor='white', linewidth=0.5)

    # Etiquetas con cobertura si > 0
    for i, (_, row) in enumerate(top.iterrows()):
        label = f"  {row['abund']:.1f}%"
        if row['cov'] > 0:
            label += f"  (cov: {row['cov']:.1f}×)"
        ax.text(row['abund'] + 0.2, i, label, va='center', fontsize=8, color='#333')

    ax.set_yticks(y)
    ax.set_yticklabels(top['Species'], fontsize=9, fontstyle='italic')
    ax.invert_yaxis()
    ax.set_xlabel('Abundancia taxonómica (%)', fontsize=10)
    ax.set_title(f'Top {top_n} especies (GTDB) — {sample_name}',
                 fontsize=12, fontweight='bold', color=BLUE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', linestyle='--', alpha=0.3)
    plt.tight_layout()
    return fig


def fig_whittaker(sample_df, sample_name):
    """Curva Whittaker (rank-abundance) — solo GTDB."""
    gtdb = get_gtdb(sample_df)
    if gtdb.empty:
        return None

    grp = gtdb.groupby('Species')['Tax_abund'].sum() \
        .sort_values(ascending=False).reset_index()
    grp['rank'] = np.arange(1, len(grp) + 1)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(grp['rank'], grp['Tax_abund'], c=BLUE, s=30, alpha=0.7, zorder=3)
    ax.plot(grp['rank'], grp['Tax_abund'], c=BLUE, alpha=0.4, linewidth=1)

    # Anotar top 5
    for _, row in grp.head(5).iterrows():
        label = row['Species'][:25] + '…' if len(row['Species']) > 27 else row['Species']
        ax.annotate(label, (row['rank'], row['Tax_abund']),
                    textcoords='offset points', xytext=(8, 4), fontsize=7.5,
                    color='#333', fontstyle='italic',
                    arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

    ax.set_yscale('log')
    ax.set_xlabel('Rango (más → menos abundante)', fontsize=10)
    ax.set_ylabel('Abundancia taxonómica (%, log)', fontsize=10)
    ax.set_title(f'Curva Rank-Abundance (Whittaker) — {sample_name}  [{len(grp)} spp.]',
                 fontsize=12, fontweight='bold', color=BLUE)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    return fig


def fig_rarefaction(sample_df, sample_name, n_steps=40, n_rep=5):
    """Curva de rarefacción simulada desde abundancias."""
    gtdb = get_gtdb(sample_df)
    if gtdb.empty or len(gtdb) < 2:
        return None

    grp = gtdb.groupby('Species')['Tax_abund'].sum().reset_index()
    grp = grp[grp['Tax_abund'] > 0]
    if len(grp) < 2:
        return None

    props = grp['Tax_abund'].values / grp['Tax_abund'].sum()
    n_total = len(grp)
    n_max = 100_000

    steps = np.unique(np.logspace(1, np.log10(n_max), n_steps).astype(int))

    means, stds = [], []
    for n in steps:
        counts = [np.sum(np.random.multinomial(n, props) > 0) for _ in range(n_rep)]
        means.append(np.mean(counts))
        stds.append(np.std(counts))

    means, stds = np.array(means), np.array(stds)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.fill_between(steps, means - stds, means + stds, alpha=0.15, color=BLUE)
    ax.plot(steps, means, color=BLUE, linewidth=2)

    ax.axhline(n_total, color='grey', linestyle=':', linewidth=1)
    ax.text(steps[-1] * 0.5, n_total * 1.03,
            f'Total: {n_total} spp.', fontsize=8, color='grey')

    # Veredicto
    plateau_pct = means[-1] / n_total * 100
    if plateau_pct > 80:
        verdict = f'Meseta alcanzada ({plateau_pct:.0f}%) — profundidad suficiente'
        vcolor = 'green'
    elif plateau_pct > 50:
        verdict = f'Curva moderada ({plateau_pct:.0f}%) — profundidad aceptable'
        vcolor = 'orange'
    else:
        verdict = f'Curva ascendente ({plateau_pct:.0f}%) — aumentar profundidad'
        vcolor = 'red'

    ax.text(0.02, 0.95, verdict, transform=ax.transAxes,
            fontsize=9, fontweight='bold', color=vcolor, va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor=vcolor, alpha=0.9))

    ax.set_xscale('log')
    ax.set_xlabel('Lecturas simuladas', fontsize=10)
    ax.set_ylabel('Especies detectadas', fontsize=10)
    ax.set_title(f'Curva de Rarefacción — {sample_name}',
                 fontsize=12, fontweight='bold', color=BLUE)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    return fig


def fig_assembly_candidates(sample_df, sample_name, min_cov=5.0):
    """
    Tabla de especies con cobertura ≥ min_cov.
    Estas son candidatas interesantes para ensamblaje dirigido.
    """
    gtdb = get_gtdb(sample_df)
    if gtdb.empty:
        return None

    grp = gtdb.groupby('Species').agg(
        abund=('Tax_abund', 'sum'),
        cov=('Eff_cov', 'mean'),
        ani=('Adj_ANI', 'mean')
    ).sort_values('cov', ascending=False).reset_index()

    candidates = grp[grp['cov'] >= min_cov].copy()
    if candidates.empty:
        return None

    n = len(candidates)
    fig, ax = plt.subplots(figsize=(14, max(3, n * 0.5 + 2)))
    ax.axis('off')

    # Header
    headers = ['Especie', 'Cob. efectiva (×)', 'ANI (%)', 'Abund. (%)', 'Viabilidad']
    col_x = [0.0, 0.45, 0.60, 0.72, 0.84]

    for j, (h, x) in enumerate(zip(headers, col_x)):
        ax.text(x, 1.0, h, transform=ax.transAxes,
                fontsize=9, fontweight='bold', color='white', va='center')
    ax.axhspan(0.94, 1.02, color=BLUE, transform=ax.transAxes, zorder=0)

    # Filas
    for i, (_, row) in enumerate(candidates.iterrows()):
        y = 0.92 - i * (0.85 / max(n, 1))
        bg = '#F0F4F7' if i % 2 == 0 else 'white'
        ax.axhspan(y - 0.02, y + 0.03, color=bg, transform=ax.transAxes, zorder=0)

        sp = row['Species'][:45] + '…' if len(row['Species']) > 47 else row['Species']
        ax.text(col_x[0], y, sp, transform=ax.transAxes,
                fontsize=8.5, fontstyle='italic', va='center')
        ax.text(col_x[1] + 0.06, y, f"{row['cov']:.1f}×",
                transform=ax.transAxes, fontsize=9, va='center', ha='center',
                fontweight='bold',
                color=GREEN if row['cov'] >= 20 else (ORANGE if row['cov'] >= 10 else '#333'))
        ax.text(col_x[2] + 0.04, y, f"{row['ani']:.1f}",
                transform=ax.transAxes, fontsize=8.5, va='center', ha='center')
        ax.text(col_x[3] + 0.04, y, f"{row['abund']:.1f}",
                transform=ax.transAxes, fontsize=8.5, va='center', ha='center')

        # Viabilidad de ensamblaje
        if row['cov'] >= 20:
            via = '●●● Excelente'
            vcolor = GREEN
        elif row['cov'] >= 10:
            via = '●●○ Buena'
            vcolor = ORANGE
        else:
            via = '●○○ Posible'
            vcolor = '#999'
        ax.text(col_x[4] + 0.06, y, via, transform=ax.transAxes,
                fontsize=8, va='center', ha='center', color=vcolor, fontweight='bold')

    ax.set_title(
        f'Candidatos a ensamblaje (cov ≥ {min_cov}×) — {sample_name}  [{n} especies]',
        fontsize=11, fontweight='bold', color=BLUE, pad=15)

    # Leyenda
    legend_text = (
        f'Criterio: cobertura efectiva ≥ {min_cov}× (estimada por Sylph).  '
        '●●● ≥20× excelente para MAG  |  ●●○ ≥10× buena  |  ●○○ ≥5× posible con más datos'
    )
    ax.text(0.5, -0.05, legend_text, transform=ax.transAxes,
            fontsize=7.5, ha='center', color='#666')

    plt.tight_layout()
    return fig


def fig_cover(tsv_path, df):
    """Portada."""
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor(BLUE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(BLUE)
    ax.axis('off')

    ax.add_patch(plt.Rectangle((0, 0), 0.015, 1, transform=ax.transAxes, color=ORANGE, zorder=5))

    ax.text(0.5, 0.80, 'Informe Metagenómico Sylph', transform=ax.transAxes,
            fontsize=26, fontweight='bold', color='white', ha='center')
    ax.text(0.5, 0.72, os.path.basename(tsv_path), transform=ax.transAxes,
            fontsize=13, color='#AACCEE', ha='center', style='italic')

    n_samples = df['Sample'].nunique()
    gtdb = get_gtdb(df)
    n_bact = gtdb['Species'].nunique()
    n_cand = len(gtdb.groupby('Species')['Eff_cov'].mean().loc[lambda x: x >= 5])

    stats = [
        ('Muestras', str(n_samples)),
        ('Spp. detectadas', str(n_bact)),
        ('Candidatos asm.', str(n_cand)),
    ]

    xs = np.linspace(0.2, 0.8, len(stats))
    for x, (label, val) in zip(xs, stats):
        rect = plt.Rectangle((x - 0.08, 0.50), 0.16, 0.12,
                              transform=ax.transAxes, color='#1A3F5C',
                              linewidth=2, edgecolor='#5B9ABF', zorder=3)
        ax.add_patch(rect)
        ax.text(x, 0.59, val, transform=ax.transAxes,
                fontsize=20, fontweight='bold', color='white', ha='center', zorder=4)
        ax.text(x, 0.51, label, transform=ax.transAxes,
                fontsize=9, color='#AACCEE', ha='center', zorder=4)

    ax.text(0.5, 0.40, 'GTDB r226 (Bacteria/Archaea) — Sylph ANI-based profiling',
            transform=ax.transAxes, fontsize=11, color='white', ha='center', fontweight='bold')
    ax.text(0.5, 0.22, f'Generado: {datetime.now().strftime("%d/%m/%Y %H:%M")}',
            transform=ax.transAxes, fontsize=10, color='#AACCEE', ha='center')
    ax.text(0.5, 0.16, 'FISABIO | EPIMOL',
            transform=ax.transAxes, fontsize=12, color='white', ha='center', fontweight='bold')

    return fig


# ─── Construcción del PDF ────────────────────────────────────────────

def build_pdf(tsv_path, out_path, min_cov=5.0):
    print(f'[INFO] Cargando: {tsv_path}')
    df = load_data(tsv_path)
    samples = sorted(df['Sample'].unique())
    gtdb = get_gtdb(df)
    print(f'[INFO] {len(samples)} muestras | {gtdb["Species"].nunique()} spp. GTDB')

    np.random.seed(42)

    with PdfPages(out_path) as pdf:
        # Portada
        fig = fig_cover(tsv_path, df)
        pdf.savefig(fig, facecolor=fig.get_facecolor())
        plt.close(fig)

        # Por muestra
        for s in samples:
            print(f'[PDF] {s}')
            sdf = df[df['Sample'] == s]

            # 1. Composición top especies
            fig = fig_composition(sdf, s)
            if fig:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            # 2. Whittaker
            fig = fig_whittaker(sdf, s)
            if fig:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            # 3. Rarefacción
            fig = fig_rarefaction(sdf, s)
            if fig:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            # 4. Candidatos a ensamblaje
            fig = fig_assembly_candidates(sdf, s, min_cov=min_cov)
            if fig:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

        # Metadata
        d = pdf.infodict()
        d['Title'] = f'Sylph Report — {os.path.basename(tsv_path)}'
        d['Author'] = 'FISABIO | EPIMOL'

    print(f'\n[OK] Informe: {out_path}')


# ─── Main ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsv', help='Fichero TSV de sylph profile')
    parser.add_argument('output', nargs='?', default=None, help='PDF de salida')
    parser.add_argument('--min-cov', type=float, default=5.0,
                        help='Cobertura mínima para candidatos a ensamblaje (default: 5×)')
    args = parser.parse_args()

    if not os.path.isfile(args.tsv):
        print(f'[ERROR] No existe: {args.tsv}')
        sys.exit(1)

    out = args.output or os.path.splitext(args.tsv)[0] + '_report.pdf'
    build_pdf(args.tsv, out, min_cov=args.min_cov)
