#!/usr/bin/env python3
"""
generate_screening_charts.py
Genera graficos de screening para EpiTax y EpiTaxMAG:
  1. FastQ Screen 100% horizontal stacked bar (estilo MultiQC/Highcharts)
  2. Sylph target organisms heatmap dual (abundancia + cobertura)
     con toggle genero/especie

Importable como modulo o ejecutable standalone.
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly

warnings.filterwarnings('ignore')

# 30 colores unicos (sin repetir) para FastQ Screen genomes
HC_PALETTE = [
    '#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
    '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1',
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#393b79',
]
GREY_MULTI = '#999999'
GREY_NOHIT = '#e6e6e6'


# ═══════════════════════════════════════════════════════════════
# FASTQ SCREEN — 100% HORIZONTAL STACKED BAR
# ═══════════════════════════════════════════════════════════════

def load_fastqscreen_raw(results_dir):
    """Load _screen.txt files directly for proper exclusive-mapping data."""
    rows = []
    for f in sorted(os.listdir(os.path.join(results_dir, '02_fastqscreen_raw')) if os.path.isdir(os.path.join(results_dir, '02_fastqscreen_raw')) else []):
        if not f.endswith('_screen.txt'):
            continue
        path = os.path.join(results_dir, '02_fastqscreen_raw', f)
        sample = f.replace('_screen.txt', '')
        try:
            with open(path) as fh:
                lines = fh.readlines()
            header_idx = None
            for i, line in enumerate(lines):
                if line.strip().startswith('Genome'):
                    header_idx = i
                    break
            if header_idx is None:
                continue
            genome_data = {}
            total_reads = 0
            for line in lines[header_idx + 1:]:
                line = line.strip()
                if not line or line.startswith('%') or line.startswith('#'):
                    break
                parts = line.split('\t')
                if len(parts) >= 12:
                    genome = parts[0]
                    reads_processed = int(parts[1])
                    if total_reads == 0:
                        total_reads = reads_processed
                    one_hit_one = int(parts[4])      # Exclusive single hit
                    multi_hit_one = int(parts[6])     # Multiple hits, one genome only
                    one_hit_multi = int(parts[8])     # One hit, multiple genomes
                    multi_hit_multi = int(parts[10])  # Multiple hits, multiple genomes

                    # Exclusive = one_hit_one_genome (truly exclusive to this genome)
                    genome_data[genome] = {
                        'exclusive': one_hit_one,
                        'multi_one': multi_hit_one,
                        'one_multi': one_hit_multi,
                        'multi_multi': multi_hit_multi,
                    }

            # Calculate categories for 100% stacked bar:
            # For each read, it goes to ONE category only:
            # 1. Exclusive to genome X (one_hit_one_genome)
            # 2. Multiple genomes (everything with multi-genome hits)
            # 3. No hits

            exclusive_total = sum(d['exclusive'] for d in genome_data.values())
            multi_genome_reads = sum(d['one_multi'] + d['multi_multi'] for d in genome_data.values())
            # multi_genome is overcounted (same read counted for each genome)
            # Use unmapped from first genome to get no_hits
            unmapped = total_reads - exclusive_total  # Simplification

            # For the stacked bar, calculate per-genome exclusive percentages
            row = {'Sample': sample, 'total_reads': total_reads}
            for genome, d in genome_data.items():
                row[genome] = d['exclusive'] / total_reads * 100 if total_reads > 0 else 0

            # "No hits" from %Hit_no_genomes line
            for line in lines:
                if line.startswith('%Hit_no_genomes'):
                    try:
                        row['No hits'] = float(line.split(':')[1].strip())
                    except:
                        pass

            # "Multiple Genomes" = 100 - exclusive_sum - no_hits
            excl_sum = sum(v for k, v in row.items() if k not in ('Sample', 'total_reads', 'No hits'))
            no_hits = row.get('No hits', 0)
            row['Multiple Genomes'] = max(0, 100 - excl_sum - no_hits)

            rows.append(row)
        except Exception:
            continue
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def build_fastqscreen_100pct(df):
    """100% Horizontal Stacked Bar Chart — MultiQC/Highcharts style."""
    if df.empty:
        return None

    samples = df['Sample'].tolist()
    genomes = [c for c in df.columns if c not in ('Sample', 'total_reads', 'No hits', 'Multiple Genomes')]

    fig = go.Figure()

    # Add genome traces with Highcharts palette cycling
    for i, genome in enumerate(genomes):
        color = HC_PALETTE[i % len(HC_PALETTE)]
        fig.add_trace(go.Bar(
            name=genome,
            y=samples,
            x=df[genome].values,
            orientation='h',
            marker_color=color,
            hovertemplate='<b>%{y}</b><br>' + genome + ': %{x:.2f}%<extra></extra>',
        ))

    # Multiple Genomes (grey)
    if 'Multiple Genomes' in df.columns:
        fig.add_trace(go.Bar(
            name='Multiple Genomes',
            y=samples,
            x=df['Multiple Genomes'].values,
            orientation='h',
            marker_color=GREY_MULTI,
            hovertemplate='<b>%{y}</b><br>Multiple Genomes: %{x:.2f}%<extra></extra>',
        ))

    # No hits (light grey)
    if 'No hits' in df.columns:
        fig.add_trace(go.Bar(
            name='No hits',
            y=samples,
            x=df['No hits'].values,
            orientation='h',
            marker_color=GREY_NOHIT,
            hovertemplate='<b>%{y}</b><br>No hits: %{x:.2f}%<extra></extra>',
        ))

    fig.update_layout(
        barmode='stack',
        title=dict(text='FastQ Screen — Distribucion de Reads Raw (%)', font=dict(size=14, color='#333')),
        xaxis=dict(title='Percentage [%]', range=[0, 100], dtick=20,
                   showgrid=True, gridcolor='#eee', gridwidth=1),
        yaxis=dict(autorange='reversed', tickfont=dict(size=9)),
        height=max(300, len(samples) * 35 + 150),
        margin=dict(l=200, r=30, t=50, b=100),
        legend=dict(font_size=8, orientation='h', yanchor='top', y=-0.15,
                    xanchor='center', x=0.5),
        font=dict(family='Helvetica, Arial, sans-serif', color='#333'),
        plot_bgcolor='white',
        bargap=0,
    )

    return fig


# ═══════════════════════════════════════════════════════════════
# SYLPH TARGET ORGANISMS — DUAL HEATMAP WITH GENUS/SPECIES TOGGLE
# ═══════════════════════════════════════════════════════════════

def load_sylph_profile(path):
    if not os.path.isfile(path):
        return pd.DataFrame()
    df = pd.read_csv(path, sep='\t')
    df['Sample'] = df['Sample_file'].apply(
        lambda s: os.path.basename(str(s)).replace('.filtered.fastq.gz', '').replace('.clean.fastq.gz', '').replace('.fastq.gz', ''))
    df['Tax_abund'] = pd.to_numeric(df['Taxonomic_abundance'], errors='coerce').fillna(0)
    df['Eff_cov'] = pd.to_numeric(df['Eff_cov'], errors='coerce').fillna(0)
    df['Contig_name'] = df['Contig_name'].astype(str)
    # Extract genus and species
    parts = df['Contig_name'].str.split(expand=True)
    df['Genus'] = parts[1].fillna('') if 1 in parts.columns else ''
    df['Species'] = (parts[1].fillna('') + ' ' + parts[2].fillna('')).str.strip() if 2 in parts.columns else df['Genus']
    # Clean MAG: prefix
    df['Genus'] = df['Genus'].str.replace('MAG:', '').str.strip()
    df['Species'] = df['Species'].str.replace('MAG:', '').str.strip()
    return df


def filter_targets(df, target_list):
    results = []
    for target in target_list:
        mask = df['Genus'].str.startswith(target, na=False)
        if target == 'Vibrio':
            mask = mask & ~df['Genus'].str.contains('Desulfovibrio', na=False)
        sub = df[mask].copy()
        if not sub.empty:
            sub['Target_genus'] = target
            results.append(sub)
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()


def build_sylph_dual_heatmap(df, target_list):
    """Build combined heatmap with buttons for Genus/Species and Abundance/Coverage."""
    if df.empty:
        return None

    filtered = filter_targets(df, target_list)
    if filtered.empty:
        return None

    # --- GENUS level aggregation ---
    genus_agg = filtered.groupby(['Sample', 'Target_genus']).agg(
        Abundance=('Tax_abund', 'sum'),
        Eff_cov=('Eff_cov', 'max')
    ).reset_index()

    pivot_g_abund = genus_agg.pivot_table(index='Target_genus', columns='Sample',
                                           values='Abundance', fill_value=0)
    pivot_g_cov = genus_agg.pivot_table(index='Target_genus', columns='Sample',
                                         values='Eff_cov', fill_value=0)

    # Sort by total abundance
    order = pivot_g_abund.sum(axis=1).sort_values(ascending=True).index
    pivot_g_abund = pivot_g_abund.loc[order]
    pivot_g_cov = pivot_g_cov.loc[order]

    # --- SPECIES level aggregation ---
    species_agg = filtered.groupby(['Sample', 'Species']).agg(
        Abundance=('Tax_abund', 'sum'),
        Eff_cov=('Eff_cov', 'max')
    ).reset_index()
    # Top 30 species
    top_sp = species_agg.groupby('Species')['Abundance'].sum().nlargest(30).index
    species_agg = species_agg[species_agg['Species'].isin(top_sp)]

    pivot_s_abund = species_agg.pivot_table(index='Species', columns='Sample',
                                             values='Abundance', fill_value=0)
    pivot_s_cov = species_agg.pivot_table(index='Species', columns='Sample',
                                           values='Eff_cov', fill_value=0)
    order_s = pivot_s_abund.sum(axis=1).sort_values(ascending=True).index
    pivot_s_abund = pivot_s_abund.loc[order_s]
    pivot_s_cov = pivot_s_cov.loc[order_s]

    # Build figure with 4 traces (2 visible at a time)
    fig = go.Figure()

    # Trace 0: Genus Abundance (VISIBLE)
    fig.add_trace(go.Heatmap(
        z=pivot_g_abund.values, x=pivot_g_abund.columns.tolist(), y=pivot_g_abund.index.tolist(),
        colorscale=[[0,'#FFFFFF'],[0.01,'#FFF3E0'],[0.1,'#FF9800'],[0.5,'#E65100'],[1,'#BF360C']],
        text=[[f'{v:.2f}%' if v>0.01 else '' for v in r] for r in pivot_g_abund.values],
        texttemplate='%{text}', textfont_size=8,
        hovertemplate='%{y} en %{x}: %{z:.3f}%<extra>Abundancia</extra>',
        colorbar=dict(title=dict(text='Abund. %', side='right'), x=1.02),
        visible=True, name='Genus - Abundancia',
    ))

    # Trace 1: Genus Coverage (HIDDEN)
    fig.add_trace(go.Heatmap(
        z=pivot_g_cov.values, x=pivot_g_cov.columns.tolist(), y=pivot_g_cov.index.tolist(),
        colorscale=[[0,'#FFFFFF'],[0.001,'#E3F2FD'],[0.01,'#42A5F5'],[0.1,'#1565C0'],[1,'#0D47A1']],
        text=[[f'{v:.2f}x' if v>0.001 else '' for v in r] for r in pivot_g_cov.values],
        texttemplate='%{text}', textfont_size=8,
        hovertemplate='%{y} en %{x}: %{z:.3f}x<extra>Cobertura</extra>',
        colorbar=dict(title=dict(text='Cob. (x)', side='right'), x=1.02),
        visible=False, name='Genus - Cobertura',
    ))

    # Trace 2: Species Abundance (HIDDEN)
    fig.add_trace(go.Heatmap(
        z=pivot_s_abund.values, x=pivot_s_abund.columns.tolist(), y=pivot_s_abund.index.tolist(),
        colorscale=[[0,'#FFFFFF'],[0.01,'#FFF3E0'],[0.1,'#FF9800'],[0.5,'#E65100'],[1,'#BF360C']],
        text=[[f'{v:.2f}%' if v>0.01 else '' for v in r] for r in pivot_s_abund.values],
        texttemplate='%{text}', textfont_size=7,
        hovertemplate='%{y} en %{x}: %{z:.3f}%<extra>Abundancia</extra>',
        colorbar=dict(title=dict(text='Abund. %', side='right'), x=1.02),
        visible=False, name='Species - Abundancia',
    ))

    # Trace 3: Species Coverage (HIDDEN)
    fig.add_trace(go.Heatmap(
        z=pivot_s_cov.values, x=pivot_s_cov.columns.tolist(), y=pivot_s_cov.index.tolist(),
        colorscale=[[0,'#FFFFFF'],[0.001,'#E3F2FD'],[0.01,'#42A5F5'],[0.1,'#1565C0'],[1,'#0D47A1']],
        text=[[f'{v:.2f}x' if v>0.001 else '' for v in r] for r in pivot_s_cov.values],
        texttemplate='%{text}', textfont_size=7,
        hovertemplate='%{y} en %{x}: %{z:.3f}x<extra>Cobertura</extra>',
        colorbar=dict(title=dict(text='Cob. (x)', side='right'), x=1.02),
        visible=False, name='Species - Cobertura',
    ))

    # Buttons: 4 options
    n_genus = len(pivot_g_abund)
    n_species = len(pivot_s_abund)
    # Use species height (largest) as fixed height so nothing gets cut
    chart_height = max(500, n_species * 22 + 180)

    fig.update_layout(
        updatemenus=[
            dict(
                type='buttons', direction='right',
                x=0.0, y=1.18, xanchor='left', yanchor='top',
                buttons=[
                    dict(label='Genero: Abundancia',
                         method='update',
                         args=[{'visible': [True, False, False, False]},
                               {'title.text': 'Alerta Temprana — Abundancia por Genero (Sylph)'}]),
                    dict(label='Genero: Cobertura',
                         method='update',
                         args=[{'visible': [False, True, False, False]},
                               {'title.text': 'Alerta Temprana — Cobertura Efectiva por Genero (Sylph)'}]),
                    dict(label='Especie: Abundancia',
                         method='update',
                         args=[{'visible': [False, False, True, False]},
                               {'title.text': 'Alerta Temprana — Abundancia por Especie (Sylph)'}]),
                    dict(label='Especie: Cobertura',
                         method='update',
                         args=[{'visible': [False, False, False, True]},
                               {'title.text': 'Alerta Temprana — Cobertura Efectiva por Especie (Sylph)'}]),
                ],
                font=dict(size=10),
                bgcolor='#f0f4f7',
                bordercolor='#2C5F8A',
            )
        ],
        title=dict(text='Alerta Temprana — Abundancia por Genero (Sylph)',
                   font=dict(size=14, color='#2C5F8A')),
        height=chart_height,
        margin=dict(l=200, r=100, t=120, b=80),
        xaxis=dict(tickangle=-40, tickfont_size=9, side='bottom'),
        yaxis=dict(tickfont_size=9),
        font=dict(family='Helvetica, Arial, sans-serif'),
    )

    return fig


# ═══════════════════════════════════════════════════════════════
# STANDALONE HTML
# ═══════════════════════════════════════════════════════════════

def build_html(fqs_chart, sylph_chart):
    plotly_js = plotly.offline.get_plotlyjs()

    def to_div(fig):
        if fig is None:
            return '<p style="color:#999;font-style:italic">Sin datos disponibles.</p>'
        return fig.to_html(full_html=False, include_plotlyjs=False)

    return f"""<!DOCTYPE html><html><head>
    <meta charset="UTF-8"><title>Screening Organismos de Inter&eacute;s</title>
    <script>{plotly_js}</script>
    <style>
    body {{ font-family: Helvetica, Arial, sans-serif; margin: 2rem; background: #f8f9fa; color: #333; }}
    .section {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 1.5rem; overflow: hidden; }}
    .section-header {{ background: #2C5F8A; color: white; padding: 0.8rem 1.5rem; font-size: 1.05rem; font-weight: 600; }}
    .section-body {{ padding: 1.5rem; }}
    .section-body > p {{ margin-bottom: 0.8rem; color: #555; font-size: 0.9rem; }}
    .note {{ font-size: 0.78rem; color: #6C757D; margin-top: 0.8rem; }}
    </style></head><body>

    <div class="section">
        <div class="section-header">FastQ Screen — Distribucion de Reads Raw</div>
        <div class="section-body">
            <p>Barras apiladas al 100% mostrando la fraccion de reads exclusivamente mapeados a cada genoma de referencia.
            Cada lectura se contabiliza una sola vez (categorias mutuamente excluyentes).
            <strong>Multiple Genomes</strong>: reads que mapean a mas de un genoma.
            <strong>No hits</strong>: reads sin mapeo a ningun genoma del panel.</p>
            {to_div(fqs_chart)}
        </div>
    </div>

    <div class="section">
        <div class="section-header">Alerta Temprana — Organismos de Inter&eacute;s (Sylph)</div>
        <div class="section-body">
            <p>Deteccion de organismos del panel de vigilancia en lecturas filtradas (Sylph).
            Usa los botones para alternar entre <strong>Genero/Especie</strong> y <strong>Abundancia/Cobertura</strong>.
            La abundancia (%) indica la fraccion del metagenoma. La cobertura efectiva (x) indica la profundidad:
            &ge;1x = deteccion fiable, &ge;30x = ensamblable.</p>
            {to_div(sylph_chart)}
            <p class="note">Lista de organismos configurable en nextflow.config (params.target_organisms).
            A nivel de genero se agregan (suma abundancia, max cobertura).
            A nivel de especie se muestran las top 30 especies detectadas.</p>
        </div>
    </div>

    </body></html>"""


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Genera graficos de screening')
    parser.add_argument('--results-dir', required=True, help='Directorio de resultados')
    parser.add_argument('--sylph', default='', help='sylph_profile_all.tsv (override)')
    parser.add_argument('--target-organisms', default='', help='Lista separada por comas')
    parser.add_argument('--output', required=True, help='Output HTML')
    args = parser.parse_args()

    rd = args.results_dir
    targets = [t.strip() for t in args.target_organisms.split(',') if t.strip()]

    # FastQ Screen from raw screen.txt files
    print('[INFO] Cargando FastQ Screen...')
    fqs_df = load_fastqscreen_raw(rd)
    fqs_chart = build_fastqscreen_100pct(fqs_df)
    print(f'[INFO] FastQ Screen: {len(fqs_df)} muestras')

    # Sylph
    sylph_path = args.sylph or f'{rd}/10_tax_sylph/sylph_profile_all.tsv'
    print(f'[INFO] Cargando Sylph: {sylph_path}')
    sylph_df = load_sylph_profile(sylph_path)
    sylph_chart = build_sylph_dual_heatmap(sylph_df, targets)
    print(f'[INFO] Sylph: {len(filter_targets(sylph_df, targets))} hits')

    html = build_html(fqs_chart, sylph_chart)
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(html)
    print(f'[OK] {args.output}')


if __name__ == '__main__':
    main()
