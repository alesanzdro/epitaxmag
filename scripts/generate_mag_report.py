#!/usr/bin/env python3
"""
generate_mag_report.py
Genera reportes HTML para EpiTaxMAG:
  - Reporte general del run (todas las muestras)
  - Reportes individuales por muestra (con branding)

Uso:
  python3 scripts/generate_mag_report.py \
      --results-dir results/260226_EPIM232 \
      --run-name 260226_EPIM232 \
      --branding branding/ \
      --output-dir results/260226_EPIM232/28_reports/
"""

import argparse, os, sys, glob, base64, warnings
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════
# PALETA Y CONSTANTES
# ═══════════════════════════════════════════════════════════════

BLUE = '#2C5F8A'; LBLUE = '#5B9ABF'; ORANGE = '#E05C2A'
GREEN = '#28A745'; RED = '#E74C3C'; GREY = '#6C757D'
PALETTE = [BLUE, ORANGE, GREEN, '#9B59B6', '#F39C12',
           '#1ABC9C', RED, '#3498DB', '#95A5A6', '#D35400']

RISK_COLORS = {'CRITICO': RED, 'ALTO': ORANGE, 'MEDIO': '#F39C12', 'BAJO': GREEN}
RISK_BADGE = {
    'CRITICO': f'<span class="badge" style="background:{RED}">CRITICO</span>',
    'ALTO':    f'<span class="badge" style="background:{ORANGE}">ALTO</span>',
    'MEDIO':   f'<span class="badge" style="background:#F39C12">MEDIO</span>',
    'BAJO':    f'<span class="badge" style="background:{GREEN}">BAJO</span>',
}
QC_BADGE = {
    'HQ': f'<span class="badge" style="background:{GREEN}">HQ</span>',
    'MQ': f'<span class="badge" style="background:{ORANGE}">MQ</span>',
    'LQ': f'<span class="badge" style="background:{GREY}">LQ</span>',
}
CONCORDANCE_BADGE = {
    'BOTH':      f'<span class="badge" style="background:{GREEN}">Reads+Contigs</span>',
    'READS':     f'<span class="badge" style="background:#F39C12">Solo Reads</span>',
    'CONTIGS':   f'<span class="badge" style="background:{LBLUE}">Solo Contigs</span>',
}

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

    for ext in ['svg', 'png', 'ico', 'jpg']:
        logo_path = os.path.join(branding_dir, f'logo.{ext}')
        if os.path.isfile(logo_path):
            mime = {'svg': 'image/svg+xml', 'png': 'image/png',
                    'ico': 'image/x-icon', 'jpg': 'image/jpeg'}[ext]
            with open(logo_path, 'rb') as f:
                logo_b64 = f'data:{mime};base64,{base64.b64encode(f.read()).decode()}'
            break

    return info, logo_b64


# ═══════════════════════════════════════════════════════════════
# LOADERS
# ═══════════════════════════════════════════════════════════════

def parse_float(val):
    if isinstance(val, (int, float)): return float(val)
    return float(str(val).replace(',', '.'))

def load_qc_csv(path):
    if not os.path.isfile(path): return pd.DataFrame()
    with open(path) as f: h = f.readline()
    sep = '\t' if '\t' in h else (';' if ';' in h else ',')
    df = pd.read_csv(path, sep=sep)
    for c in df.columns:
        if c != 'Sample':
            df[c] = df[c].apply(lambda v: parse_float(v) if pd.notna(v) else 0)
    return df

def load_run_stats(input_dir):
    """Calcula estadisticas de la carrera incluyendo unclassified."""
    import gzip
    stats = {'samples': [], 'total_reads': 0, 'classified_reads': 0,
             'unclassified_reads': 0, 'pct_unclassified': 0}
    if not input_dir or not os.path.isdir(input_dir):
        return stats
    for f in sorted(glob.glob(f'{input_dir}/*.fastq.gz')):
        name = os.path.basename(f).replace('.fastq.gz','')
        try:
            n = 0
            with gzip.open(f, 'rt') as fh:
                for line in fh:
                    if line.startswith('@'): n += 1
            stats['total_reads'] += n
            if 'unclassified' in name.lower():
                stats['unclassified_reads'] = n
            else:
                stats['classified_reads'] += n
                stats['samples'].append({'name': name, 'reads': n})
        except: pass
    if stats['total_reads'] > 0:
        stats['pct_unclassified'] = round(stats['unclassified_reads'] / stats['total_reads'] * 100, 1)
    return stats

def load_checkm2(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/23_binqc_checkm2/*/quality_report.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                comp = float(r['Completeness']); cont = float(r['Contamination'])
                qc = 'HQ' if comp >= 90 and cont <= 5 else ('MQ' if comp >= 50 and cont <= 10 else 'LQ')
                rows.append({'Sample': sample, 'Bin': r['Name'], 'Completeness': comp,
                             'Contamination': cont, 'Quality': qc})
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_taxonomy(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/24_taxonomy_gtdbtk/*/*_taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                cls = str(r.iloc[1]) if len(r) > 1 else ''
                levels = {}
                for part in cls.split(';'):
                    p = part.strip()
                    for px, lv in [('d__','domain'),('p__','phylum'),('c__','class'),
                                   ('o__','order'),('f__','family'),('g__','genus'),('s__','species')]:
                        if p.startswith(px): levels[lv] = p[len(px):] or ''
                org = f"{levels.get('genus','')} {levels.get('species','')}".strip() or 'Unclassified'
                rows.append({'Sample': sample, 'Bin': str(r.iloc[0]).strip(),
                             'Taxonomy': cls, 'Organism': org,
                             'Genus': levels.get('genus',''), 'Species': levels.get('species','')})
        except: pass
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
                        'Sample': sample, 'Bin': str(r.get('Name','')).strip(),
                        'Contig': str(r.get('Contig id','')).strip(),
                        'Gene': str(r.get('Element symbol', r.get('Gene symbol',''))).strip(),
                        'Gene_desc': str(r.get('Element name', r.get('Sequence name',''))).strip(),
                        'Scope': str(r.get('Scope','')).strip(),
                        'Type': str(r.get('Type', r.get('Element type',''))).strip(),
                        'Subtype': str(r.get('Subtype', r.get('Element subtype',''))).strip(),
                        'Class': str(r.get('Class','')).strip(),
                        'Subclass': str(r.get('Subclass','')).strip(),
                        'Method': str(r.get('Method','')).strip(),
                        'Identity': float(r.get('% Identity to reference', r.get('% Identity to reference sequence', 0))),
                        'Coverage': float(r.get('% Coverage of reference', r.get('% Coverage of reference sequence', 0))),
                        'Accession': str(r.get('Closest reference accession','')).strip(),
                        'Closest_ref': str(r.get('Closest reference name','')).strip(),
                        'HMM_desc': str(r.get('HMM description','')).strip(),
                        'Start': str(r.get('Start','')).strip(),
                        'Stop': str(r.get('Stop','')).strip(),
                        'Strand': str(r.get('Strand','')).strip(),
                    })
            except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_kma(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/11_amr_kma/*.res')):
        sample = os.path.basename(f).replace('.res','')
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                template = str(r.iloc[0]).strip()
                if not template or template.startswith('#'): continue
                gene_root = template.split('_')[0] if '_' in template else template
                rows.append({
                    'Sample': sample, 'Template': template, 'Gene_root': gene_root,
                    'KMA_Identity': float(r.get('Template_Identity', 0)),
                    'KMA_Coverage': float(r.get('Template_Coverage', 0)),
                    'KMA_Depth': float(r.get('Depth', 0)),
                })
        except: pass
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
                rows.append({'Sample': sample, 'Contig': contig, 'Plasmid_score': score})
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_integration(results_dir):
    dfs = []
    for f in sorted(glob.glob(f'{results_dir}/27_amr_pathogen_integration/*_amr_pathogen.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty: dfs.append(df)
        except: pass
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


# ═══════════════════════════════════════════════════════════════
# ANALYSIS
# ═══════════════════════════════════════════════════════════════

def compute_survival(qc_raw, qc_filt):
    if qc_raw.empty or qc_filt.empty: return pd.DataFrame()
    raw = qc_raw[['Sample','total_reads','mean_length']].copy()
    raw.columns = ['Sample','raw_reads','raw_mean_len']
    raw['raw_gb'] = raw['raw_reads'] * raw['raw_mean_len'] / 1e9
    filt = qc_filt[['Sample','total_reads','mean_length']].copy()
    filt.columns = ['Sample','clean_reads','clean_mean_len']
    filt['clean_gb'] = filt['clean_reads'] * filt['clean_mean_len'] / 1e9
    df = pd.merge(raw, filt, on='Sample', how='outer').fillna(0)
    df['pct_reads'] = np.where(df['raw_reads']>0, df['clean_reads']/df['raw_reads']*100, 0).round(1)
    df['pct_gb'] = np.where(df['raw_gb']>0, df['clean_gb']/df['raw_gb']*100, 0).round(1)
    return df.sort_values('Sample').reset_index(drop=True)

def compute_concordance(kma_df, amr_df):
    """Cruzar KMA reads vs AMRFinderPlus contigs por gene root."""
    if kma_df.empty and amr_df.empty: return pd.DataFrame()
    rows = []
    samples = set()
    if not kma_df.empty: samples.update(kma_df['Sample'].unique())
    if not amr_df.empty: samples.update(amr_df['Sample'].unique())

    for sample in sorted(samples):
        kma_genes = set()
        kma_info = {}
        if not kma_df.empty:
            sk = kma_df[kma_df['Sample']==sample]
            for _, r in sk.iterrows():
                root = r['Gene_root']
                kma_genes.add(root)
                if root not in kma_info or r['KMA_Depth'] > kma_info[root]['KMA_Depth']:
                    kma_info[root] = {'KMA_Depth': r['KMA_Depth'],
                                      'KMA_Identity': r['KMA_Identity'],
                                      'KMA_Coverage': r['KMA_Coverage']}

        mag_genes = set()
        mag_info = {}
        if not amr_df.empty:
            sa = amr_df[amr_df['Sample']==sample]
            for _, r in sa.iterrows():
                gene = r['Gene']
                mag_genes.add(gene)
                if gene not in mag_info:
                    mag_info[gene] = {'MAG_Identity': r['Identity'], 'MAG_Coverage': r['Coverage'],
                                      'Bin': r['Bin'], 'Class': r['Class']}

        all_genes = kma_genes | mag_genes
        for gene in sorted(all_genes):
            in_kma = gene in kma_genes
            in_mag = gene in mag_genes
            source = 'BOTH' if in_kma and in_mag else ('READS' if in_kma else 'CONTIGS')
            row = {'Sample': sample, 'Gene': gene, 'Source': source}
            if in_kma:
                row.update(kma_info.get(gene, {}))
            if in_mag:
                row.update(mag_info.get(gene, {}))
            rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame()

def confidence_score(identity, coverage, method, kma_concordance, mag_quality, plasmid_score):
    """Score de confianza compuesto (0-100)."""
    score = 0
    score += min(identity, 100) * 0.3
    score += min(coverage, 100) * 0.2
    method_scores = {'EXACT': 20, 'ALLELE': 18, 'BLAST': 15, 'PARTIAL': 10,
                     'HMM': 12, 'INTERNAL_STOP': 8}
    score += method_scores.get(method, 10)
    if kma_concordance: score += 15
    qc_scores = {'HQ': 15, 'MQ': 10, 'LQ': 5}
    score += qc_scores.get(mag_quality, 0)
    return min(round(score), 100)


# ═══════════════════════════════════════════════════════════════
# PLOTLY CHARTS
# ═══════════════════════════════════════════════════════════════

def html_pipeline_dag():
    """Genera diagrama del pipeline estilo nf-core con versiones."""
    return """<div class="dag-container"><div class="dag-flow">
    <div class="dag-row">
        <div class="dag-box qc">NanoPlot 1.46.2</div>
        <div class="dag-box qc">FastQC 0.12.1</div>
        <div class="dag-box qc">CHECK_SQCORE</div>
        <div class="dag-box qc">FastQ Screen 0.16.0</div>
    </div>
    <div class="dag-label">QC + Screening (reads crudos)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box trim">Porechop ABI 0.5.1</div>
        <div class="dag-arrow">&#x25B6;</div>
        <div class="dag-box trim">Chopper 0.12.0</div>
        <div class="dag-arrow">&#x25B6;</div>
        <div class="dag-box trim">Host Removal (minimap2 2.30)</div>
    </div>
    <div class="dag-label">Trimming + Filtrado + Decontaminacion</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box qc">NanoPlot 1.46.2</div>
        <div class="dag-box qc">FastQC 0.12.1</div>
        <div class="dag-box qc">CHECK_SQCORE</div>
    </div>
    <div class="dag-label">QC reads filtrados</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box tax">Kraken2 2.17.1</div>
        <div class="dag-box tax">Bracken 3.1</div>
        <div class="dag-box tax">Kaiju 1.10.1</div>
        <div class="dag-box tax">Sylph 0.9.0</div>
        <div class="dag-box ann">KMA 1.6.8</div>
    </div>
    <div class="dag-label">EpiTax: Taxonomia + AMR reads (ResFinder)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box rep">MultiQC 1.33</div>
        <div class="dag-box rep">EpiTax HTML Report</div>
    </div>
    <div class="dag-label">Reportes fase TAX</div>
    <div class="dag-arrow">&#x25BC; reads filtrados</div>

    <div class="dag-row">
        <div class="dag-box asm">MetaFlye 2.9.6</div>
        <div class="dag-arrow">&#x25B6;</div>
        <div class="dag-box asm">Medaka 2.2.1</div>
        <div class="dag-arrow">&#x25B6;</div>
        <div class="dag-box qc">QUAST 5.3.0</div>
    </div>
    <div class="dag-label">Ensamblaje + Polishing</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box asm">minimap2 2.30 + samtools 1.21</div>
    </div>
    <div class="dag-label">Mapeo cobertura</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box bin">MetaBAT2 2.18</div>
        <div class="dag-box bin">MaxBin2 2.2.7</div>
        <div class="dag-box bin">SemiBin2 2.2.1</div>
    </div>
    <div class="dag-label">Binning (3 algoritmos)</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box bin">DAS Tool 1.1.7</div>
    </div>
    <div class="dag-label">Refinamiento de bins</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box qc">CheckM2 1.1.0</div>
        <div class="dag-box tax">GTDB-Tk 2.7.0 (r232)</div>
        <div class="dag-box ann">AMRFinderPlus 4.2.7</div>
        <div class="dag-box ann">geNomad 1.12.0</div>
    </div>
    <div class="dag-label">QC + Taxonomia + AMR contigs + Plasmidos</div>
    <div class="dag-arrow">&#x25BC;</div>

    <div class="dag-row">
        <div class="dag-box rep">Integracion AMR-Patogeno-Plasmido</div>
        <div class="dag-arrow">&#x25B6;</div>
        <div class="dag-box rep">EpiTaxMAG Report</div>
    </div>
    <div class="dag-label">Reporte final con evaluacion de riesgo</div>
    </div></div>"""

def fig_to_html(fig):
    if fig is None: return '<p class="no-data">Sin datos disponibles.</p>'
    return fig.to_html(full_html=False, include_plotlyjs=False)

def chart_retention(surv):
    if surv.empty: return None
    fig = make_subplots(rows=1, cols=2, subplot_titles=['Reads','Gigabases'], horizontal_spacing=0.12)
    s = surv['Sample'].values
    fig.add_trace(go.Bar(name='Raw', x=s, y=surv['raw_reads'], marker_color=BLUE, opacity=0.7,
                         text=[f'{v:,.0f}' for v in surv['raw_reads']], textposition='outside', textfont_size=8), row=1, col=1)
    fig.add_trace(go.Bar(name='Clean', x=s, y=surv['clean_reads'], marker_color=GREEN, opacity=0.85,
                         text=[f'{v:,.0f}' for v in surv['clean_reads']], textposition='outside', textfont_size=8), row=1, col=1)
    fig.add_trace(go.Bar(name='Raw Gb', x=s, y=surv['raw_gb'], marker_color=BLUE, opacity=0.7, showlegend=False,
                         text=[f'{v:.2f}' for v in surv['raw_gb']], textposition='outside', textfont_size=8), row=1, col=2)
    fig.add_trace(go.Bar(name='Clean Gb', x=s, y=surv['clean_gb'], marker_color=GREEN, opacity=0.85, showlegend=False,
                         text=[f'{v:.2f}' for v in surv['clean_gb']], textposition='outside', textfont_size=8), row=1, col=2)
    fig.update_layout(barmode='group', height=400, margin=dict(l=60,r=30,t=60,b=100),
                      legend=dict(orientation='h', yanchor='bottom', y=1.04, xanchor='center', x=0.5),
                      font=dict(family='Arial, sans-serif'))
    fig.update_xaxes(tickangle=-40, tickfont_size=8)
    return fig

def chart_mag_quality(checkm2_df):
    if checkm2_df.empty: return None
    counts = checkm2_df.groupby(['Sample','Quality']).size().unstack(fill_value=0)
    for q in ['HQ','MQ','LQ']:
        if q not in counts.columns: counts[q] = 0
    fig = go.Figure()
    colors = {'HQ': GREEN, 'MQ': ORANGE, 'LQ': GREY}
    for q in ['HQ','MQ','LQ']:
        fig.add_trace(go.Bar(name=q, x=counts.index, y=counts[q], marker_color=colors[q],
                             text=counts[q], textposition='inside'))
    fig.update_layout(barmode='stack', height=350, title=dict(text='Calidad de MAGs por muestra', font_color=BLUE),
                      xaxis_tickangle=-40, margin=dict(l=60,r=30,t=60,b=100),
                      legend=dict(orientation='h', yanchor='bottom', y=1.04, xanchor='center', x=0.5),
                      font=dict(family='Arial, sans-serif', size=10))
    return fig

def chart_amr_heatmap(integ_df):
    if integ_df.empty: return None
    pivot = integ_df.groupby(['Organism','AMR_class']).size().unstack(fill_value=0)
    if pivot.empty: return None
    top_org = pivot.sum(axis=1).nlargest(15).index
    top_cls = pivot.sum(axis=0).nlargest(10).index
    pivot = pivot.loc[pivot.index.isin(top_org), pivot.columns.isin(top_cls)]
    if pivot.empty: return None
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=[[0,'#FFFFFF'],[0.01,'#F0F4F7'],[0.1,'#5B9ABF'],[0.5,'#2C5F8A'],[1,'#1A3F5C']],
        text=[[str(v) if v > 0 else '' for v in row] for row in pivot.values],
        texttemplate='%{text}', textfont_size=9,
        colorbar=dict(title=dict(text='N genes', side='right'))))
    fig.update_layout(title=dict(text='AMR: Organismos vs Clases de Resistencia', font_color=BLUE, font_size=14),
                      height=max(350, len(pivot)*28+150), margin=dict(l=250,r=80,t=60,b=150),
                      xaxis=dict(tickangle=-45, tickfont_size=8), yaxis=dict(tickfont_size=9),
                      font=dict(family='Arial, sans-serif'))
    return fig

def chart_risk_summary(integ_df):
    if integ_df.empty: return None
    risk_counts = integ_df['Risk_level'].value_counts()
    colors_map = [RISK_COLORS.get(r, GREY) for r in risk_counts.index]
    fig = go.Figure(go.Bar(x=risk_counts.index, y=risk_counts.values, marker_color=colors_map,
                           text=risk_counts.values, textposition='outside'))
    fig.update_layout(title=dict(text='Distribucion de Riesgo AMR', font_color=BLUE, font_size=14),
                      height=300, margin=dict(l=60,r=30,t=60,b=60),
                      yaxis_title='N genes', font=dict(family='Arial, sans-serif'))
    return fig


# ═══════════════════════════════════════════════════════════════
# HTML TABLES
# ═══════════════════════════════════════════════════════════════

def html_survival_table(surv):
    if surv.empty: return '<p class="no-data">Sin datos de QC.</p>'
    rows = ''
    for _, r in surv.iterrows():
        pct_cls = 'good' if r['pct_reads'] >= 70 else ('warn' if r['pct_reads'] >= 50 else 'bad')
        gb_cls = 'highlight' if r['clean_gb'] >= 1.0 else ''
        rows += (f'<tr><td class="sample" data-sample="{r["Sample"]}">{r["Sample"]}</td>'
                 f'<td>{r["raw_reads"]:,.0f}</td><td>{r["raw_gb"]:.2f}</td>'
                 f'<td>{r["clean_reads"]:,.0f}</td><td class="{gb_cls}">{r["clean_gb"]:.2f}</td>'
                 f'<td class="{pct_cls}">{r["pct_reads"]:.1f}%</td></tr>\n')
    tot = surv[['raw_reads','raw_gb','clean_reads','clean_gb']].sum()
    pct = (tot['clean_reads']/tot['raw_reads']*100) if tot['raw_reads']>0 else 0
    rows += (f'<tr class="total-row"><td><strong>TOTAL</strong></td>'
             f'<td><strong>{tot["raw_reads"]:,.0f}</strong></td><td><strong>{tot["raw_gb"]:.2f}</strong></td>'
             f'<td><strong>{tot["clean_reads"]:,.0f}</strong></td><td><strong>{tot["clean_gb"]:.2f}</strong></td>'
             f'<td><strong>{pct:.1f}%</strong></td></tr>')
    return (f'<table class="data-table sortable filterable"><thead><tr><th>Muestra</th><th>Raw Reads</th><th>Raw Gb</th>'
            f'<th>Clean Reads</th><th>Clean Gb</th><th>% Retencion</th></tr></thead><tbody>{rows}</tbody></table>')

def html_mag_catalog(checkm2_df, tax_df, sample=None):
    if checkm2_df.empty: return '<p class="no-data">Sin MAGs recuperados.</p>'
    df = checkm2_df.copy()
    if sample: df = df[df['Sample']==sample]
    if not tax_df.empty:
        df = df.merge(tax_df[['Sample','Bin','Organism','Taxonomy']], on=['Sample','Bin'], how='left')
    else:
        df['Organism'] = ''; df['Taxonomy'] = ''
    df = df.fillna('')
    rows = ''
    for _, r in df.iterrows():
        sample_td = f'<td class="sample" data-sample="{r.get("Sample","")}">{r.get("Sample","")}</td>' if not sample else ''
        rows += (f'<tr>{sample_td}<td>{r["Bin"]}</td>'
                 f'<td><em>{r.get("Organism","")}</em></td>'
                 f'<td>{r["Completeness"]:.1f}%</td><td>{r["Contamination"]:.1f}%</td>'
                 f'<td>{QC_BADGE.get(r["Quality"],"")}</td></tr>\n')
    sample_th = '<th>Muestra</th>' if not sample else ''
    return (f'<table class="data-table sortable filterable"><thead><tr>{sample_th}'
            f'<th>Bin</th><th>Organismo</th>'
            f'<th>Completeness</th><th>Contamination</th><th>Calidad</th></tr></thead>'
            f'<tbody>{rows}</tbody></table>')

def _safe_str(val):
    """Return empty string for nan/NA/None, otherwise str."""
    if val is None: return ''
    if isinstance(val, float) and np.isnan(val): return ''
    s = str(val).strip()
    return '' if s.lower() in ('nan', 'na', 'none', '') else s

def html_amr_detail(integ_df, amr_mags_df, sample=None):
    """Enhanced AMR detail table with descriptions, accession links, and closest reference."""
    if integ_df.empty: return '<p class="no-data">Sin genes AMR detectados en MAGs.</p>'
    df = integ_df.copy()
    if sample: df = df[df['Sample']==sample]
    if df.empty: return '<p class="no-data">Sin genes AMR en esta muestra.</p>'

    # Build lookup from amr_mags for extra columns (Accession, Closest_ref, HMM_desc)
    amr_lookup = {}
    if not amr_mags_df.empty:
        for _, ar in amr_mags_df.iterrows():
            key = (str(ar.get('Sample','')), str(ar.get('Bin','')),
                   str(ar.get('Gene','')), str(ar.get('Contig','')))
            amr_lookup[key] = {
                'Accession': _safe_str(ar.get('Accession','')),
                'Closest_ref': _safe_str(ar.get('Closest_ref','')),
                'HMM_desc': _safe_str(ar.get('HMM_desc','')),
            }

    rows = ''
    for _, r in df.iterrows():
        conf = confidence_score(r.get('Identity_pct',0), r.get('Coverage_pct',0),
                                r.get('Method',''), r.get('Location','')=='PLASMID',
                                r.get('MAG_quality',''), r.get('Plasmid_score',0))
        conf_cls = 'good' if conf >= 70 else ('warn' if conf >= 40 else 'bad')
        loc = r.get('Location','')
        if loc == 'PLASMID':
            loc_badge = f'<span class="badge" style="background:{RED}">PLASMIDO</span>'
        else:
            loc_badge = f'<span class="badge" style="background:{LBLUE}">CROMOSOMA</span>'
        s_name = str(r.get("Sample",""))
        sample_td = f'<td class="sample" data-sample="{s_name}">{s_name}</td>' if not sample else ''

        # Lookup extra AMR info
        gene_sym = str(r.get('Gene',''))
        mag_name = str(r.get('MAG', r.get('Bin','')))
        contig_name = str(r.get('Contig',''))
        lookup_key = (s_name, mag_name, gene_sym, contig_name)
        extra = amr_lookup.get(lookup_key, {})
        accession = extra.get('Accession','')
        closest_ref = extra.get('Closest_ref','')
        # Fallback: also try without contig (less specific)
        if not accession and not closest_ref:
            for k, v in amr_lookup.items():
                if k[0] == s_name and k[2] == gene_sym:
                    accession = v.get('Accession','')
                    closest_ref = v.get('Closest_ref','')
                    break

        # Description: use Closest_ref (has the specific variant), fallback to Gene_description
        description = closest_ref if closest_ref else str(r.get('Gene_description',''))

        # Build links
        links_parts = []
        if accession:
            ncbi_url = f'https://www.ncbi.nlm.nih.gov/protein/{accession}'
            links_parts.append(f'<a href="{ncbi_url}" target="_blank" title="NCBI Protein">{accession}</a>')
        if gene_sym:
            card_url = f'https://card.mcmaster.ca/ontology/query?query={gene_sym}'
            links_parts.append(f'<a href="{card_url}" target="_blank" title="CARD Database">CARD</a>')
        links_html = ' | '.join(links_parts) if links_parts else '-'

        method_val = str(r.get('Method',''))

        rows += (f'<tr>{sample_td}'
                 f'<td><em>{r.get("Organism","")}</em></td>'
                 f'<td><strong>{gene_sym}</strong></td>'
                 f'<td class="desc-cell">{description}</td>'
                 f'<td>{r.get("AMR_class","")}</td>'
                 f'<td>{loc_badge}</td>'
                 f'<td>{r.get("Identity_pct",0):.1f}%</td>'
                 f'<td>{r.get("Coverage_pct",0):.1f}%</td>'
                 f'<td>{method_val}</td>'
                 f'<td class="{conf_cls}"><strong>{conf}</strong></td>'
                 f'<td>{RISK_BADGE.get(r.get("Risk_level",""),"")}</td>'
                 f'<td class="links-cell">{links_html}</td></tr>\n')

    sample_th = '<th>Muestra</th>' if not sample else ''
    return (f'<table class="data-table sortable filterable"><thead><tr>'
            f'{sample_th}'
            f'<th>Organismo</th><th>Gen</th><th>Descripcion</th><th>Clase AMR</th><th>Ubicacion</th>'
            f'<th>% Identidad</th><th>% Cobertura</th><th>Metodo</th>'
            f'<th>Confianza</th><th>Riesgo</th><th>Links</th></tr></thead>'
            f'<tbody>{rows}</tbody></table>'
            f'<p class="note">Confianza: score compuesto (0-100) basado en identidad, cobertura, '
            f'metodo de deteccion, concordancia KMA-MAG y calidad del MAG. '
            f'Descripcion: referencia mas cercana de AMRFinderPlus (variante especifica). '
            f'Links: NCBI Protein (accession) y CARD (base de datos de resistencia).</p>')

def html_concordance_table(conc_df, sample=None):
    if conc_df.empty: return '<p class="no-data">Sin datos de concordancia KMA-MAG.</p>'
    df = conc_df.copy()
    if sample: df = df[df['Sample']==sample]
    if df.empty: return '<p class="no-data">Sin datos de concordancia para esta muestra.</p>'
    rows = ''
    for _, r in df.iterrows():
        badge = CONCORDANCE_BADGE.get(r.get('Source',''), '')
        kma_d = f"{r.get('KMA_Depth',0):.1f}x" if pd.notna(r.get('KMA_Depth')) and r.get('KMA_Depth',0)>0 else '-'
        kma_i = f"{r.get('KMA_Identity',0):.1f}%" if pd.notna(r.get('KMA_Identity')) and r.get('KMA_Identity',0)>0 else '-'
        mag_i = f"{r.get('MAG_Identity',0):.1f}%" if pd.notna(r.get('MAG_Identity')) and r.get('MAG_Identity',0)>0 else '-'
        mag_b = str(r.get('Bin','')) if pd.notna(r.get('Bin')) else '-'
        s_name = str(r.get('Sample',''))
        sample_td = f'<td class="sample" data-sample="{s_name}">{s_name}</td>' if not sample else ''
        rows += (f'<tr>{sample_td}'
                 f'<td><strong>{r["Gene"]}</strong></td>'
                 f'<td>{r.get("Class","")}</td>'
                 f'<td>{badge}</td>'
                 f'<td>{kma_d}</td><td>{kma_i}</td>'
                 f'<td>{mag_i}</td><td>{mag_b}</td></tr>\n')
    sample_th = '<th>Muestra</th>' if not sample else ''
    return (f'<table class="data-table sortable filterable"><thead><tr>'
            f'{sample_th}'
            f'<th>Gen</th><th>Clase</th><th>Deteccion</th>'
            f'<th>KMA Depth</th><th>KMA Identity</th>'
            f'<th>MAG Identity</th><th>MAG (Bin)</th></tr></thead>'
            f'<tbody>{rows}</tbody></table>'
            f'<p class="note">Reads+Contigs = maxima confianza. Solo Reads = posible reservorio no ensamblado. '
            f'Solo Contigs = profundidad de reads insuficiente para KMA.</p>')

def html_plasmid_amr(integ_df, sample=None):
    if integ_df.empty: return '<p class="no-data">Sin datos de plasmidos.</p>'
    df = integ_df[integ_df['Location']=='PLASMID'].copy()
    if sample: df = df[df['Sample']==sample]
    if df.empty: return '<p class="no-data">Sin genes AMR en plasmidos.</p>'
    rows = ''
    for _, r in df.iterrows():
        s_name = str(r.get('Sample',''))
        sample_td = f'<td class="sample" data-sample="{s_name}">{s_name}</td>' if not sample else ''
        rows += (f'<tr>{sample_td}'
                 f'<td><em>{r.get("Organism","")}</em></td>'
                 f'<td><strong>{r.get("Gene","")}</strong></td>'
                 f'<td>{r.get("AMR_class","")}</td>'
                 f'<td>{r.get("Contig","")}</td>'
                 f'<td><strong>{r.get("Plasmid_score",0):.2f}</strong></td>'
                 f'<td>{RISK_BADGE.get(r.get("Risk_level",""),"")}</td></tr>\n')
    sample_th = '<th>Muestra</th>' if not sample else ''
    return (f'<table class="data-table sortable filterable"><thead><tr>'
            f'{sample_th}'
            f'<th>Organismo</th><th>Gen</th><th>Clase AMR</th>'
            f'<th>Contig</th><th>Score Plasmido</th><th>Riesgo</th></tr></thead>'
            f'<tbody>{rows}</tbody></table>'
            f'<p class="note">Score plasmido (geNomad): &ge;0.9 seguro, 0.7-0.9 probable, &lt;0.7 dudoso. '
            f'Genes AMR en plasmidos conjugativos representan riesgo de diseminacion horizontal.</p>')


# ═══════════════════════════════════════════════════════════════
# CSS - separate <style> block only
# ═══════════════════════════════════════════════════════════════

HTML_CSS = """<style>
:root { --blue: #2C5F8A; --orange: #E05C2A; --green: #28A745; --red: #E74C3C; --grey: #6C757D; --bg: #F8F9FA; }
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: 'Segoe UI', Tahoma, Arial, sans-serif; background: var(--bg); color: #333; line-height: 1.6; }
.header { background: linear-gradient(135deg, var(--blue) 0%, #1a3f5c 100%); color: white; padding: 1.5rem 2rem; border-bottom: 4px solid var(--orange); display: flex; align-items: center; gap: 1.5rem; }
.header-logo { height: 60px; }
.header-text h1 { font-size: 1.5rem; margin-bottom: 0.2rem; }
.header-text .subtitle { opacity: 0.8; font-size: 0.85rem; }
.header-text .org-info { opacity: 0.6; font-size: 0.75rem; }
.container { max-width: 1200px; margin: 0 auto; padding: 1.5rem; }
.metrics-row { display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 1.5rem; }
.metric-card { flex: 1; min-width: 150px; background: white; border-radius: 8px; padding: 1rem; box-shadow: 0 2px 6px rgba(0,0,0,0.08); text-align: center; }
.metric-card .value { font-size: 1.8rem; font-weight: 700; color: var(--blue); }
.metric-card .label { font-size: 0.8rem; color: var(--grey); }
.section { background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 1.5rem; overflow: hidden; }
.section-header { background: var(--blue); color: white; padding: 0.8rem 1.5rem; font-size: 1.05rem; font-weight: 600; }
.section-body { padding: 1.5rem; }
.section-body > p { margin-bottom: 0.8rem; color: #555; font-size: 0.9rem; }
.data-table { width: 100%; border-collapse: collapse; font-size: 0.82rem; margin: 1rem 0; }
.data-table th { background: var(--blue); color: white; padding: 0.5rem 0.6rem; text-align: center; font-weight: 600; font-size: 0.78rem; cursor: pointer; position: relative; user-select: none; }
.data-table th:hover { background: #1a3f5c; }
.data-table td { padding: 0.4rem 0.6rem; text-align: center; border-bottom: 1px solid #eee; }
.data-table tr:nth-child(even) { background: #F0F4F7; }
.data-table tr:hover { background: #E8EFF5; }
.data-table .sample { text-align: left; font-weight: 600; font-family: monospace; font-size: 0.78rem; }
.data-table .good { color: var(--green); font-weight: 600; }
.data-table .warn { color: var(--orange); font-weight: 600; }
.data-table .bad { color: var(--red); font-weight: 600; }
.data-table .highlight { background: #C8E6C9 !important; font-weight: 700; }
.data-table .total-row { background: #E3EBF3 !important; font-size: 0.85rem; }
.data-table .desc-cell { text-align: left; font-size: 0.76rem; max-width: 220px; word-wrap: break-word; }
.data-table .links-cell { font-size: 0.74rem; white-space: nowrap; }
.data-table .links-cell a { color: var(--blue); text-decoration: none; }
.data-table .links-cell a:hover { text-decoration: underline; }
.badge { display: inline-block; padding: 0.15rem 0.5rem; border-radius: 10px; font-size: 0.7rem; font-weight: 700; color: white; }
.note { font-size: 0.78rem; color: var(--grey); margin-top: 0.8rem; }
.no-data { color: #999; font-style: italic; padding: 1rem 0; }
.alert-box { background: #FFEBEE; border-left: 4px solid var(--red); padding: 1rem; border-radius: 4px; margin: 1rem 0; }
.alert-box p { margin-bottom: 0.5rem; }
.conclusion { background: #E8F5E9; border-left: 4px solid var(--green); padding: 1rem; border-radius: 4px; }
.footer { text-align: center; padding: 1.5rem; color: var(--grey); font-size: 0.75rem; }
.sample-filter { background: #f0f4f7; padding: 0.8rem 1.2rem; border-radius: 6px; margin-bottom: 1rem; }
.sample-filter label { margin-right: 1rem; font-size: 0.85rem; cursor: pointer; }
.sample-filter input[type=checkbox] { margin-right: 0.3rem; }
.search-box { margin-bottom: 1rem; }
.search-box input { width: 100%; padding: 0.5rem 0.8rem; border: 1px solid #ccc; border-radius: 4px; font-size: 0.85rem; }
.search-box input:focus { outline: none; border-color: var(--blue); box-shadow: 0 0 0 2px rgba(44,95,138,0.2); }
.sort-asc::after { content: ' \\25B2' !important; opacity: 1 !important; }
.sort-desc::after { content: ' \\25BC' !important; opacity: 1 !important; }
.data-table.sortable th::after { content: ' \\2195'; font-size: 0.65rem; opacity: 0.5; }
.dag-container { padding: 1rem; overflow-x: auto; }
.dag-flow { display: flex; flex-direction: column; gap: 0.3rem; align-items: center; font-size: 0.75rem; }
.dag-row { display: flex; gap: 0.5rem; align-items: center; flex-wrap: wrap; justify-content: center; }
.dag-box { padding: 0.4rem 0.7rem; border-radius: 4px; color: white; font-weight: 600; white-space: nowrap; text-align: center; }
.dag-box.qc { background: #5B9ABF; }
.dag-box.trim { background: #2C5F8A; }
.dag-box.tax { background: #28A745; }
.dag-box.asm { background: #E05C2A; }
.dag-box.bin { background: #9B59B6; }
.dag-box.ann { background: #E74C3C; }
.dag-box.rep { background: #1ABC9C; }
.dag-arrow { color: #aaa; font-size: 1.2rem; }
.dag-label { font-size: 0.65rem; color: #999; margin-top: -0.2rem; }
@media print { .section { break-inside: avoid; } body { background: white; } .sample-filter, .search-box { display: none; } }
</style>"""


# ═══════════════════════════════════════════════════════════════
# JS - separate <script> block, placed before </body>
# ═══════════════════════════════════════════════════════════════

HTML_JS = """<script>
(function() {
    'use strict';

    // ── Table sorting: click any <th> in a table with class "sortable" ──
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
            // Determine sort direction
            var asc = !th.classList.contains('sort-asc');
            // Clear all sort indicators in this table
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
                // Try numeric comparison (strip %, commas, x suffix)
                var na = parseFloat(va.replace(/[,%x]/g, '').replace(/[^0-9.\\-]/g, ''));
                var nb = parseFloat(vb.replace(/[,%x]/g, '').replace(/[^0-9.\\-]/g, ''));
                if (!isNaN(na) && !isNaN(nb)) {
                    return asc ? na - nb : nb - na;
                }
                // Alphabetic fallback
                return asc ? va.localeCompare(vb) : vb.localeCompare(va);
            });
            for (var j = 0; j < rows.length; j++) {
                tbody.appendChild(rows[j]);
            }
        });
    }

    // ── Search box: filters ALL tables with class "filterable" ──
    function initSearchFilter() {
        var input = document.getElementById('searchAmr');
        if (!input) return;
        input.addEventListener('input', function() {
            var filter = this.value.toLowerCase();
            var tables = document.querySelectorAll('table.filterable');
            for (var t = 0; t < tables.length; t++) {
                var trs = tables[t].querySelectorAll('tbody tr');
                for (var i = 0; i < trs.length; i++) {
                    var text = trs[i].textContent.toLowerCase();
                    trs[i].style.display = (filter === '' || text.indexOf(filter) !== -1) ? '' : 'none';
                }
            }
        });
    }

    // ── Sample filter checkboxes ──
    // Toggles visibility of <tr> rows where any <td> has data-sample matching
    // Note: this affects TABLES only, NOT Plotly charts (charts are rendered as SVG/canvas)
    window.toggleSample = function(sampleName, checked) {
        var tds = document.querySelectorAll('td[data-sample="' + sampleName + '"]');
        for (var i = 0; i < tds.length; i++) {
            var row = tds[i].closest('tr');
            if (row) {
                row.style.display = checked ? '' : 'none';
            }
        }
    };

    // ── Initialize on DOM ready ──
    document.addEventListener('DOMContentLoaded', function() {
        initSortableTables();
        initSearchFilter();
    });
})();
</script>"""


# ═══════════════════════════════════════════════════════════════
# HTML BUILDERS
# ═══════════════════════════════════════════════════════════════

def build_sample_checkboxes(surv):
    if surv.empty: return ''
    items = []
    for s in sorted(surv['Sample'].unique()):
        if 'unclassified' in s.lower(): continue
        items.append(f'<label><input type="checkbox" checked '
                     f"onchange=\"toggleSample('{s}', this.checked)\"> {s}</label>")
    return ' '.join(items)

def build_header(org_info, logo_b64, title, subtitle, date_str):
    logo_html = f'<img src="{logo_b64}" class="header-logo" alt="Logo">' if logo_b64 else ''
    org_line = f'{org_info.get("name","")} | {org_info.get("department","")}'
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
        EpiTaxMAG Pipeline v1.0 | {org_info.get('name','')} | {date_str}
    </div>"""

def build_general_report(data, org_info, logo_b64, run_name):
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')
    plotly_js = plotly.offline.get_plotlyjs()
    surv = data['survival']
    checkm2 = data['checkm2']
    integ = data['integration']
    conc = data['concordance']
    amr_mags = data['amr_mags']
    run_stats = data.get('run_stats')

    # Metrics cards
    n_samples = len(surv) if not surv.empty else 0
    total_raw_gb = surv['raw_gb'].sum() if not surv.empty else 0
    total_clean_gb = surv['clean_gb'].sum() if not surv.empty else 0
    mean_ret = surv['pct_reads'].mean() if not surv.empty else 0
    std_ret = surv['pct_reads'].std() if not surv.empty and len(surv)>1 else 0
    n_mags = len(checkm2) if not checkm2.empty else 0
    n_hq = len(checkm2[checkm2['Quality']=='HQ']) if not checkm2.empty else 0
    n_amr = len(integ) if not integ.empty else 0
    n_plasmid = len(integ[integ['Location']=='PLASMID']) if not integ.empty else 0
    n_critical = len(integ[integ['Risk_level']=='CRITICO']) if not integ.empty else 0

    # Unclassified stats
    pct_unclass = run_stats['pct_unclassified'] if run_stats else 0
    total_run_reads = run_stats['total_reads'] if run_stats else 0
    unclass_reads = run_stats['unclassified_reads'] if run_stats else 0
    unclass_color = RED if pct_unclass >= 25 else (ORANGE if pct_unclass >= 15 else GREEN)
    unclass_warn = ''
    if pct_unclass >= 25:
        unclass_warn = f'<span style="color:{RED}; font-weight:700"> (ALERTA: &ge;25%)</span>'
    elif pct_unclass >= 15:
        unclass_warn = f'<span style="color:{ORANGE}; font-weight:700"> (Atencion: &ge;15%)</span>'

    metrics_html = f"""<div class="metrics-row">
        <div class="metric-card"><div class="value">{n_samples}</div><div class="label">Muestras</div></div>
        <div class="metric-card"><div class="value">{total_raw_gb:.1f}</div><div class="label">Gb Raw Total</div></div>
        <div class="metric-card"><div class="value">{total_clean_gb:.1f}</div><div class="label">Gb Clean Total</div></div>
        <div class="metric-card"><div class="value">{mean_ret:.1f}%</div><div class="label">Retencion ({chr(177)}{std_ret:.1f}%)</div></div>
        <div class="metric-card"><div class="value" style="color:{unclass_color}">{pct_unclass}%</div><div class="label">Unclassified</div></div>
        <div class="metric-card"><div class="value">{n_mags}</div><div class="label">MAGs ({n_hq} HQ)</div></div>
        <div class="metric-card"><div class="value">{n_amr}</div><div class="label">Genes AMR</div></div>
        <div class="metric-card"><div class="value" style="color:{RED if n_critical>0 else GREEN}">{n_critical}</div><div class="label">Riesgo Critico</div></div>
    </div>"""

    # Alerts
    alerts = ''
    if n_critical > 0:
        alert_items = ''
        for _, r in integ[integ['Risk_level']=='CRITICO'].iterrows():
            alert_items += f'<p><strong>{r.get("Gene","")}</strong> ({r.get("AMR_class","")}) en <em>{r.get("Organism","")}</em> - plasmido score {r.get("Plasmid_score",0):.2f} [{r.get("Sample","")}]</p>'
        alerts = f'<div class="alert-box"><p><strong>ALERTA: {n_critical} genes AMR criticos en plasmidos (transferibles)</strong></p>{alert_items}</div>'

    # Charts
    retention_chart = fig_to_html(chart_retention(surv))
    quality_chart = fig_to_html(chart_mag_quality(checkm2))
    heatmap_chart = fig_to_html(chart_amr_heatmap(integ))
    risk_chart = fig_to_html(chart_risk_summary(integ))

    header = build_header(org_info, logo_b64, f'EpiTaxMAG Report - {run_name}',
                          f'Reporte General | {n_samples} muestras', date_str)
    footer = build_footer(org_info, date_str)

    return f"""<!DOCTYPE html><html lang="es"><head>
    <meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EpiTaxMAG Report - {run_name}</title>
    <script type="text/javascript">{plotly_js}</script>
    {HTML_CSS}
    </head><body>
    {header}
    <div class="container">
        {metrics_html}
        {alerts}

        <div class="section">
            <div class="section-header">Pipeline EpiTaxMAG - Flujo de Trabajo</div>
            <div class="section-body">
                <p>Diagrama completo del pipeline con versiones de software utilizadas.</p>
                {html_pipeline_dag()}
            </div>
        </div>

        <div class="sample-filter" id="sampleFilter">
            <strong>Filtrar muestras:</strong>
            {build_sample_checkboxes(surv)}
            <p class="note" style="margin-top:0.4rem">Nota: los filtros de muestra solo afectan a las tablas, no a los graficos Plotly.</p>
        </div>

        <div class="section">
            <div class="section-header">0. Estadisticas de la Carrera de Secuenciacion</div>
            <div class="section-body">
                <p>Metricas globales de la carrera Nanopore. El porcentaje de lecturas <strong>unclassified</strong>
                (sin barcode asignado) indica el aprovechamiento del demultiplexado.
                Un valor &ge;25% es preocupante y puede indicar problemas con el kit de barcodes o la calidad de la libreria.</p>
                <table class="data-table">
                    <thead><tr><th>Metrica</th><th>Valor</th><th>Evaluacion</th></tr></thead>
                    <tbody>
                        <tr><td>Total reads secuenciados</td><td><strong>{total_run_reads:,}</strong></td><td>-</td></tr>
                        <tr><td>Reads clasificados (con barcode)</td><td>{total_run_reads - unclass_reads:,}</td><td>-</td></tr>
                        <tr><td>Reads unclassified (sin barcode)</td><td>{unclass_reads:,}</td>
                            <td style="color:{unclass_color}; font-weight:700">{pct_unclass}%{unclass_warn}</td></tr>
                        <tr><td>Muestras demultiplexadas</td><td>{n_samples}</td><td>-</td></tr>
                        <tr><td>Gb raw total (muestras)</td><td>{total_raw_gb:.2f}</td><td>-</td></tr>
                        <tr><td>Gb clean total (tras filtrado)</td><td>{total_clean_gb:.2f}</td><td>-</td></tr>
                        <tr><td>Retencion media (reads)</td><td>{mean_ret:.1f}% ({chr(177)}{std_ret:.1f}%)</td>
                            <td style="color:{'var(--green)' if mean_ret>=70 else 'var(--orange)'}">
                            {'Buena' if mean_ret>=70 else 'Revisar filtros'}</td></tr>
                    </tbody>
                </table>
            </div>
        </div>

        <div class="section">
            <div class="section-header">1. Retencion de Lecturas por Muestra</div>
            <div class="section-body">
                <p>Comparativa raw vs clean tras Porechop ABI (trimming) + Chopper (filtrado Q/longitud). Celdas verdes: &ge;1 Gb limpio.</p>
                {html_survival_table(surv)}
                {retention_chart}
            </div>
        </div>

        <div class="section">
            <div class="section-header">2. Recuperacion de MAGs</div>
            <div class="section-body">
                <p>MAGs recuperados por MetaBAT2 + MaxBin2 + SemiBin2, refinados por DAS Tool, evaluados por CheckM2 y clasificados por GTDB-Tk (r232).</p>
                {html_mag_catalog(checkm2, data['taxonomy'])}
                {quality_chart}
            </div>
        </div>

        <div class="search-box">
            <input type="text" id="searchAmr" placeholder="Buscar gen, organismo, clase AMR...">
        </div>

        <div class="section">
            <div class="section-header">3. Resistencia Antimicrobiana (AMR) en MAGs</div>
            <div class="section-body">
                <p>Genes AMR detectados por AMRFinderPlus en MAGs. Ubicacion cromosoma/plasmido determinada por geNomad. Confianza: score compuesto (identidad + cobertura + metodo + concordancia + calidad MAG).</p>
                {html_amr_detail(integ, amr_mags)}
                {heatmap_chart}
            </div>
        </div>

        <div class="section">
            <div class="section-header">4. Concordancia KMA Reads vs MAG Contigs</div>
            <div class="section-body">
                <p>Comparacion de genes AMR detectados en reads (KMA/ResFinder, fase TAX) vs contigs ensamblados (AMRFinderPlus, fase MAG). Deteccion en ambos = maxima confianza. Solo en reads = posible reservorio no ensamblado.</p>
                {html_concordance_table(conc)}
            </div>
        </div>

        <div class="section">
            <div class="section-header">5. Genes AMR en Plasmidos</div>
            <div class="section-body">
                <p>Genes AMR localizados en contigs clasificados como plasmidicos por geNomad. Representan riesgo de transferencia horizontal entre organismos.</p>
                {html_plasmid_amr(integ)}
            </div>
        </div>

        <div class="section">
            <div class="section-header">6. Evaluacion de Riesgo</div>
            <div class="section-body">
                {risk_chart}
                <div class="note" style="margin-top:1rem">
                    <strong>Niveles de riesgo:</strong><br>
                    CRITICO = resistencia clave (beta-lactam, quinolona, colistina) en plasmido<br>
                    ALTO = cualquier AMR en plasmido (transferible)<br>
                    MEDIO = resistencia clave en cromosoma<br>
                    BAJO = AMR no critica en cromosoma
                </div>
            </div>
        </div>
    </div>
    {footer}
    {HTML_JS}
    </body></html>"""


def build_sample_report(sample, data, org_info, logo_b64, run_name):
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')
    plotly_js = plotly.offline.get_plotlyjs()
    surv = data['survival']
    checkm2 = data['checkm2']
    integ = data['integration']
    conc = data['concordance']
    amr_mags = data['amr_mags']
    run_stats = data.get('run_stats')

    # Sample-level data
    s_surv = surv[surv['Sample']==sample].iloc[0] if not surv.empty and sample in surv['Sample'].values else None
    s_checkm2 = checkm2[checkm2['Sample']==sample] if not checkm2.empty else pd.DataFrame()
    s_integ = integ[integ['Sample']==sample] if not integ.empty else pd.DataFrame()
    s_conc = conc[conc['Sample']==sample] if not conc.empty else pd.DataFrame()
    s_amr_mags = amr_mags[amr_mags['Sample']==sample] if not amr_mags.empty else pd.DataFrame()

    # Run context
    n_samples = len(surv) if not surv.empty else 0
    total_clean = surv['clean_gb'].sum() if not surv.empty else 0
    mean_ret = surv['pct_reads'].mean() if not surv.empty else 0
    pct_unclass = run_stats['pct_unclassified'] if run_stats else 0

    # Sample metrics
    raw_reads = f'{s_surv["raw_reads"]:,.0f}' if s_surv is not None else '-'
    raw_gb = f'{s_surv["raw_gb"]:.2f}' if s_surv is not None else '-'
    clean_reads = f'{s_surv["clean_reads"]:,.0f}' if s_surv is not None else '-'
    clean_gb = f'{s_surv["clean_gb"]:.2f}' if s_surv is not None else '-'
    pct_ret = f'{s_surv["pct_reads"]:.1f}%' if s_surv is not None else '-'
    n_bins = len(s_checkm2)
    n_hq = len(s_checkm2[s_checkm2['Quality']=='HQ']) if not s_checkm2.empty else 0
    n_amr = len(s_integ)
    n_plas = len(s_integ[s_integ['Location']=='PLASMID']) if not s_integ.empty else 0
    n_crit = len(s_integ[s_integ['Risk_level']=='CRITICO']) if not s_integ.empty else 0

    metrics_html = f"""
    <div class="metrics-row">
        <div class="metric-card"><div class="value">{raw_gb}</div><div class="label">Gb Raw</div></div>
        <div class="metric-card"><div class="value">{clean_gb}</div><div class="label">Gb Clean</div></div>
        <div class="metric-card"><div class="value">{pct_ret}</div><div class="label">Retencion</div></div>
        <div class="metric-card"><div class="value">{n_bins}</div><div class="label">MAGs ({n_hq} HQ)</div></div>
        <div class="metric-card"><div class="value">{n_amr}</div><div class="label">Genes AMR</div></div>
        <div class="metric-card"><div class="value">{n_plas}</div><div class="label">En Plasmido</div></div>
        <div class="metric-card"><div class="value" style="color:{RED if n_crit>0 else GREEN}">{n_crit}</div><div class="label">Criticos</div></div>
    </div>"""

    # Alert
    alerts = ''
    if n_crit > 0:
        items = ''
        for _, r in s_integ[s_integ['Risk_level']=='CRITICO'].iterrows():
            items += f'<p><strong>{r.get("Gene","")}</strong> ({r.get("AMR_class","")}) en <em>{r.get("Organism","")}</em> - plasmido score {r.get("Plasmid_score",0):.2f}</p>'
        alerts = f'<div class="alert-box"><p><strong>ALERTA: {n_crit} genes AMR criticos</strong></p>{items}</div>'

    header = build_header(org_info, logo_b64, f'EpiTaxMAG - {sample}',
                          f'Reporte Individual | Run: {run_name} ({n_samples} muestras, {total_clean:.1f} Gb clean, {mean_ret:.1f}% retencion, {pct_unclass}% unclassified)',
                          date_str)
    footer = build_footer(org_info, date_str)

    return f"""<!DOCTYPE html><html lang="es"><head>
    <meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EpiTaxMAG - {sample}</title>
    <script type="text/javascript">{plotly_js}</script>
    {HTML_CSS}
    </head><body>
    {header}
    <div class="container">
        {metrics_html}
        {alerts}

        <div class="section">
            <div class="section-header">1. Catalogo de MAGs</div>
            <div class="section-body">
                <p>MAGs recuperados, clasificados taxonomicamente por GTDB-Tk (r232) y evaluados por CheckM2.</p>
                {html_mag_catalog(s_checkm2, data['taxonomy'], sample)}
            </div>
        </div>

        <div class="search-box">
            <input type="text" id="searchAmr" placeholder="Buscar gen, organismo, clase AMR...">
        </div>

        <div class="section">
            <div class="section-header">2. Resistencia Antimicrobiana Detallada</div>
            <div class="section-body">
                <p>Genes AMR detectados por AMRFinderPlus. Ubicacion (cromosoma/plasmido) por geNomad. Score de confianza compuesto.</p>
                {html_amr_detail(s_integ, s_amr_mags, sample)}
            </div>
        </div>

        <div class="section">
            <div class="section-header">3. Concordancia KMA Reads vs MAG Contigs</div>
            <div class="section-body">
                <p>Deteccion cruzada: reads (KMA) vs contigs (AMRFinderPlus). Concordancia = maxima confianza. Solo reads = posible reservorio no ensamblado.</p>
                {html_concordance_table(s_conc, sample)}
            </div>
        </div>

        <div class="section">
            <div class="section-header">4. Genes AMR en Plasmidos</div>
            <div class="section-body">
                <p>AMR en elementos moviles. Score geNomad &ge;0.9 = plasmido confirmado.</p>
                {html_plasmid_amr(s_integ, sample)}
            </div>
        </div>
    </div>
    {footer}
    {HTML_JS}
    </body></html>"""


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Genera reportes HTML EpiTaxMAG')
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--run-name', required=True)
    parser.add_argument('--input-dir', default=None, help='Dir con fastq.gz originales (para stats de unclassified)')
    parser.add_argument('--branding', default='branding/')
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()

    rd = args.results_dir
    print(f'[INFO] Generando reportes EpiTaxMAG para {args.run_name}')

    org_info, logo_b64 = load_branding(args.branding)
    print(f'[INFO] Branding: {org_info.get("name","")}')

    # Run stats (unclassified)
    run_stats = None
    if args.input_dir:
        print(f'[INFO] Calculando estadisticas de carrera desde {args.input_dir}...')
        run_stats = load_run_stats(args.input_dir)
        print(f'[INFO] Carrera: {run_stats["total_reads"]:,} reads totales, '
              f'{run_stats["unclassified_reads"]:,} unclassified ({run_stats["pct_unclassified"]}%)')

    # Load all data
    qc_raw = load_qc_csv(f'{rd}/01_qc_raw/qc_raw_metrics.csv')
    qc_filt = load_qc_csv(f'{rd}/06_qc_filtered/qc_filtered_metrics.csv')
    survival = compute_survival(qc_raw, qc_filt)
    checkm2 = load_checkm2(rd)
    taxonomy = load_taxonomy(rd)
    amr_mags = load_amr_mags(rd)
    kma = load_kma(rd)
    plasmids = load_plasmids(rd)
    integration = load_integration(rd)
    concordance = compute_concordance(kma, amr_mags)

    print(f'[INFO] Datos: {len(survival)} muestras, {len(checkm2)} MAGs, '
          f'{len(integration)} AMR integrados, {len(kma)} KMA hits')

    data = {'survival': survival, 'checkm2': checkm2, 'taxonomy': taxonomy,
            'amr_mags': amr_mags, 'kma': kma, 'integration': integration,
            'concordance': concordance, 'run_stats': run_stats}

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f'{args.output_dir}/per_sample', exist_ok=True)

    # General report
    general_html = build_general_report(data, org_info, logo_b64, args.run_name)
    general_path = f'{args.output_dir}/{args.run_name}_epitaxmag_report.html'
    with open(general_path, 'w', encoding='utf-8') as f:
        f.write(general_html)
    print(f'[OK] Reporte general: {general_path}')

    # Per-sample reports
    samples = sorted(set())
    if not survival.empty: samples = sorted(survival['Sample'].unique())
    elif not checkm2.empty: samples = sorted(checkm2['Sample'].unique())

    for sample in samples:
        if 'unclassified' in sample.lower():
            continue
        sample_html = build_sample_report(sample, data, org_info, logo_b64, args.run_name)
        sample_path = f'{args.output_dir}/per_sample/{sample}_report.html'
        with open(sample_path, 'w', encoding='utf-8') as f:
            f.write(sample_html)
        print(f'[OK] Reporte muestra: {sample_path}')

    # Export combined TSV for external analysis
    if not integration.empty:
        tsv_path = f'{args.output_dir}/{args.run_name}_amr_pathogen_matrix.tsv'
        integration.to_csv(tsv_path, sep='\t', index=False)
        print(f'[OK] TSV exportado: {tsv_path}')

    print(f'\n[DONE] {len(samples)+1} reportes generados en {args.output_dir}/')


if __name__ == '__main__':
    main()
