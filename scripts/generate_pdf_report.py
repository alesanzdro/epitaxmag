#!/usr/bin/env python3
"""
generate_pdf_report.py
Genera reporte PDF profesional completo para EpiTaxMAG.

Secciones:
  A. Resumen ejecutivo
  B. Estadisticas de carrera
  C. Retencion de lecturas
  D. Catalogo de MAGs
  E. AMR detallado con confianza
  F. Concordancia KMA-MAG
  G. Virulencia (VFDB)
  H. Plasmidos y movilidad
  I. Integrones
  J. Evaluacion de riesgo
  K. Desglose por muestra
  L. Software y versiones

Uso:
  python3 scripts/generate_pdf_report.py \
      --results-dir results/260226_EPIM232 \
      --run-name 260226_EPIM232 \
      --input-dir /path/to/raw/fastqs \
      --branding branding/ \
      --output results/260226_EPIM232/32_reports/260226_EPIM232_report.pdf
"""

import argparse, os, glob, warnings, gzip, textwrap
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io

from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm, mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, PageBreak, Image, HRFlowable,
                                 KeepTogether, ListFlowable, ListItem)
from reportlab.platypus.tableofcontents import TableOfContents

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════
# COLORS
# ═══════════════════════════════════════════════════════════════
C_BLUE    = colors.HexColor('#2C5F8A')
C_LBLUE   = colors.HexColor('#5B9ABF')
C_ORANGE  = colors.HexColor('#E05C2A')
C_GREEN   = colors.HexColor('#28A745')
C_RED     = colors.HexColor('#E74C3C')
C_GREY    = colors.HexColor('#6C757D')
C_LGREY   = colors.HexColor('#F0F4F7')
C_WHITE   = colors.white

# ═══════════════════════════════════════════════════════════════
# STYLES
# ═══════════════════════════════════════════════════════════════
def get_styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle('DocTitle', parent=styles['Title'], fontSize=22,
                               textColor=C_BLUE, spaceAfter=4, alignment=TA_CENTER))
    styles.add(ParagraphStyle('DocSub', parent=styles['Normal'], fontSize=10,
                               textColor=C_GREY, alignment=TA_CENTER, spaceAfter=12))
    styles.add(ParagraphStyle('SectionH', parent=styles['Heading1'], fontSize=14,
                               textColor=C_BLUE, spaceBefore=16, spaceAfter=8,
                               borderWidth=1, borderColor=C_BLUE, borderPadding=4))
    styles.add(ParagraphStyle('SubH', parent=styles['Heading2'], fontSize=11,
                               textColor=C_LBLUE, spaceBefore=10, spaceAfter=4))
    styles.add(ParagraphStyle('Body', parent=styles['Normal'], fontSize=9,
                               textColor=colors.black, alignment=TA_JUSTIFY,
                               spaceAfter=6, leading=12))
    styles.add(ParagraphStyle('Small', parent=styles['Normal'], fontSize=7.5,
                               textColor=C_GREY, spaceAfter=4))
    styles.add(ParagraphStyle('Alert', parent=styles['Normal'], fontSize=9,
                               textColor=C_RED, fontName='Helvetica-Bold', spaceAfter=6))
    styles.add(ParagraphStyle('CellH', parent=styles['Normal'], fontSize=7.5,
                               textColor=C_WHITE, fontName='Helvetica-Bold',
                               alignment=TA_CENTER))
    styles.add(ParagraphStyle('Cell', parent=styles['Normal'], fontSize=7,
                               alignment=TA_CENTER))
    styles.add(ParagraphStyle('CellL', parent=styles['Normal'], fontSize=7,
                               alignment=TA_LEFT))
    return styles

# ═══════════════════════════════════════════════════════════════
# BRANDING
# ═══════════════════════════════════════════════════════════════
def load_branding(branding_dir):
    info = {'name': 'EpiTaxMAG', 'department': '', 'address': '',
            'contact': '', 'phone': '', 'web': ''}
    logo_path = None
    info_path = os.path.join(branding_dir, 'org.info')
    if os.path.isfile(info_path):
        with open(info_path) as f:
            for line in f:
                if '=' in line:
                    k, v = line.strip().split('=', 1)
                    info[k.strip()] = v.strip()
    for ext in ['png', 'jpg', 'svg']:
        lp = os.path.join(branding_dir, f'logo.{ext}')
        if os.path.isfile(lp):
            logo_path = lp
            break
    return info, logo_path

# ═══════════════════════════════════════════════════════════════
# DATA LOADERS (reuse same logic)
# ═══════════════════════════════════════════════════════════════
def pf(val):
    if isinstance(val, (int, float)): return float(val)
    return float(str(val).replace(',', '.'))

def load_qc(path):
    if not os.path.isfile(path): return pd.DataFrame()
    with open(path) as f: h = f.readline()
    sep = '\t' if '\t' in h else (';' if ';' in h else ',')
    df = pd.read_csv(path, sep=sep)
    for c in df.columns:
        if c != 'Sample': df[c] = df[c].apply(lambda v: pf(v) if pd.notna(v) else 0)
    return df

def load_run_stats(input_dir):
    stats = {'total_reads': 0, 'classified_reads': 0, 'unclassified_reads': 0, 'pct_unclassified': 0}
    if not input_dir or not os.path.isdir(input_dir): return stats
    for f in sorted(glob.glob(f'{input_dir}/*.fastq.gz')):
        name = os.path.basename(f).replace('.fastq.gz', '')
        try:
            n = sum(1 for line in gzip.open(f, 'rt') if line.startswith('@'))
            stats['total_reads'] += n
            if 'unclassified' in name.lower(): stats['unclassified_reads'] = n
            else: stats['classified_reads'] += n
        except: pass
    if stats['total_reads'] > 0:
        stats['pct_unclassified'] = round(stats['unclassified_reads'] / stats['total_reads'] * 100, 1)
    return stats

def load_all_data(rd):
    d = {}
    # QC
    d['qc_raw'] = load_qc(f'{rd}/01_qc_raw/qc_raw_metrics.csv')
    d['qc_filt'] = load_qc(f'{rd}/06_qc_filtered/qc_filtered_metrics.csv')
    # Survival
    if not d['qc_raw'].empty and not d['qc_filt'].empty:
        raw = d['qc_raw'][['Sample','total_reads','mean_length']].copy()
        raw.columns = ['Sample','raw_reads','raw_ml']
        raw['raw_gb'] = raw['raw_reads'] * raw['raw_ml'] / 1e9
        filt = d['qc_filt'][['Sample','total_reads','mean_length']].copy()
        filt.columns = ['Sample','clean_reads','clean_ml']
        filt['clean_gb'] = filt['clean_reads'] * filt['clean_ml'] / 1e9
        surv = pd.merge(raw, filt, on='Sample', how='outer').fillna(0)
        surv['pct'] = np.where(surv['raw_reads']>0, surv['clean_reads']/surv['raw_reads']*100, 0).round(1)
        d['survival'] = surv.sort_values('Sample').reset_index(drop=True)
    else:
        d['survival'] = pd.DataFrame()
    # CheckM2
    rows = []
    for f in sorted(glob.glob(f'{rd}/23_binqc_checkm2/*/quality_report.tsv')):
        s = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                c = float(r['Completeness']); ct = float(r['Contamination'])
                q = 'HQ' if c>=90 and ct<=5 else ('MQ' if c>=50 and ct<=10 else 'LQ')
                rows.append({'Sample':s,'Bin':r['Name'],'Comp':c,'Cont':ct,'Q':q})
        except: pass
    d['checkm2'] = pd.DataFrame(rows) if rows else pd.DataFrame()
    # Taxonomy
    rows = []
    for f in sorted(glob.glob(f'{rd}/24_taxonomy_gtdbtk/*/*_taxonomy.tsv')):
        s = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                cls = str(r.iloc[1]) if len(r)>1 else ''
                g = ''; sp = ''
                for p in cls.split(';'):
                    if p.strip().startswith('g__'): g = p.strip()[3:]
                    if p.strip().startswith('s__'): sp = p.strip()[3:]
                org = f"{g} {sp}".strip() or 'Unclassified'
                rows.append({'Sample':s,'Bin':str(r.iloc[0]).strip(),'Tax':cls,'Org':org})
        except: pass
    d['taxonomy'] = pd.DataFrame(rows) if rows else pd.DataFrame()
    # Integration
    dfs = []
    for f in sorted(glob.glob(f'{rd}/27_amr_pathogen_integration/*_amr_pathogen.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty: dfs.append(df)
        except: pass
    d['integration'] = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    # KMA
    rows = []
    for f in sorted(glob.glob(f'{rd}/11_amr_kma/*.res')):
        s = os.path.basename(f).replace('.res','')
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                t = str(r.iloc[0]).strip()
                if not t or t.startswith('#'): continue
                rows.append({'Sample':s,'Template':t,'Gene':t.split('_')[0],
                             'Depth':float(r.get('Depth',0)),
                             'Identity':float(r.get('Template_Identity',0)),
                             'Coverage':float(r.get('Template_Coverage',0))})
        except: pass
    d['kma'] = pd.DataFrame(rows) if rows else pd.DataFrame()
    # VFDB
    rows = []
    for f in sorted(glob.glob(f'{rd}/28_abricate/*/*.vfdb.tsv')):
        s = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            for _, r in df.iterrows():
                rows.append({'Sample':s,'Bin':os.path.basename(f).replace('.vfdb.tsv',''),
                             'Gene':str(r.get('GENE','')),'Product':str(r.get('PRODUCT','')),
                             'Identity':float(r.get('%IDENTITY',0)),'Coverage':float(r.get('%COVERAGE',0))})
        except: pass
    d['vfdb'] = pd.DataFrame(rows) if rows else pd.DataFrame()
    # Plasmids
    rows = []
    for f in sorted(glob.glob(f'{rd}/26_genomad/*/*_plasmid_summary.tsv')):
        s = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            rows.append({'Sample':s,'N_plasmids':len(df)})
        except: pass
    d['plasmid_counts'] = pd.DataFrame(rows) if rows else pd.DataFrame()
    return d

# ═══════════════════════════════════════════════════════════════
# TABLE HELPER
# ═══════════════════════════════════════════════════════════════
def make_table(headers, data_rows, col_widths=None):
    """Create a styled ReportLab table."""
    st = get_styles()
    tdata = [[Paragraph(h, st['CellH']) for h in headers]]
    for row in data_rows:
        tdata.append([Paragraph(str(c), st['Cell']) if i > 0 else Paragraph(str(c), st['CellL'])
                       for i, c in enumerate(row)])
    t = Table(tdata, colWidths=col_widths, repeatRows=1)
    style = [
        ('BACKGROUND', (0,0), (-1,0), C_BLUE),
        ('TEXTCOLOR', (0,0), (-1,0), C_WHITE),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('ALIGN', (0,0), (0,-1), 'LEFT'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('FONTSIZE', (0,0), (-1,-1), 7),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [C_LGREY, C_WHITE]),
        ('GRID', (0,0), (-1,-1), 0.3, colors.HexColor('#CCCCCC')),
        ('TOPPADDING', (0,0), (-1,-1), 2),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2),
        ('LEFTPADDING', (0,0), (-1,-1), 3),
        ('RIGHTPADDING', (0,0), (-1,-1), 3),
    ]
    t.setStyle(TableStyle(style))
    return t

# ═══════════════════════════════════════════════════════════════
# PAGE TEMPLATE (header/footer with branding)
# ═══════════════════════════════════════════════════════════════
class BrandedDocTemplate(SimpleDocTemplate):
    def __init__(self, filename, org_info, logo_path, run_name, **kwargs):
        self.org_info = org_info
        self.logo_path = logo_path
        self.run_name = run_name
        self.date_str = datetime.now().strftime('%d/%m/%Y')
        super().__init__(filename, **kwargs)

    def afterPage(self):
        c = self.canv
        w, h = self.pagesize
        # Footer
        c.setFont('Helvetica', 7)
        c.setFillColor(colors.HexColor('#999999'))
        footer = f"{self.org_info.get('name','')} | {self.org_info.get('department','')} | {self.org_info.get('address','')}"
        c.drawString(1.5*cm, 1*cm, footer)
        c.drawRightString(w - 1.5*cm, 1*cm, f"EpiTaxMAG | {self.run_name} | {self.date_str} | Pag. {c.getPageNumber()}")
        # Header line
        c.setStrokeColor(C_BLUE)
        c.setLineWidth(1.5)
        c.line(1.5*cm, h - 1.3*cm, w - 1.5*cm, h - 1.3*cm)
        # Logo
        if self.logo_path:
            try:
                c.drawImage(self.logo_path, 1.5*cm, h - 1.2*cm, width=1.5*cm,
                           height=1*cm, preserveAspectRatio=True, mask='auto')
            except: pass

# ═══════════════════════════════════════════════════════════════
# BUILD PDF
# ═══════════════════════════════════════════════════════════════
def build_pdf(data, run_stats, org_info, logo_path, run_name, output_path):
    st = get_styles()
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')
    surv = data['survival']
    checkm2 = data['checkm2']
    tax = data['taxonomy']
    integ = data['integration']
    kma = data['kma']
    vfdb = data['vfdb']

    doc = BrandedDocTemplate(output_path, org_info, logo_path, run_name,
                              pagesize=A4, rightMargin=1.5*cm, leftMargin=1.5*cm,
                              topMargin=2*cm, bottomMargin=1.8*cm,
                              title=f'EpiTaxMAG Report - {run_name}',
                              author=org_info.get('name', ''))

    story = []
    pw = A4[0] - 3*cm  # page width usable

    # ── COVER ──────────────────────────────────────────────────
    story.append(Spacer(1, 3*cm))
    if logo_path:
        story.append(Image(logo_path, width=4*cm, height=3*cm))
        story.append(Spacer(1, 0.5*cm))
    story.append(Paragraph('EpiTaxMAG', st['DocTitle']))
    story.append(Paragraph('Reporte de Vigilancia Metagenomica', ParagraphStyle(
        'SubTitle', parent=st['Normal'], fontSize=14, textColor=C_LBLUE, alignment=TA_CENTER, spaceAfter=6)))
    story.append(HRFlowable(width='60%', thickness=2, color=C_ORANGE, spaceAfter=12))
    story.append(Paragraph(f'Run: <b>{run_name}</b>', ParagraphStyle(
        'RunInfo', parent=st['Normal'], fontSize=11, alignment=TA_CENTER, spaceAfter=4)))
    story.append(Paragraph(f'Fecha: {date_str}', st['DocSub']))
    story.append(Paragraph(f'{org_info.get("name","")} | {org_info.get("department","")}', st['DocSub']))
    story.append(Paragraph(f'{org_info.get("address","")}', st['DocSub']))

    # Quick stats on cover
    n_samples = len(surv) if not surv.empty else 0
    n_mags = len(checkm2) if not checkm2.empty else 0
    n_amr = len(integ) if not integ.empty else 0
    n_crit = len(integ[integ['Risk_level']=='CRITICO']) if not integ.empty else 0

    story.append(Spacer(1, 1*cm))
    cover_data = [
        ['Muestras', str(n_samples)],
        ['Total Gb (raw)', f"{surv['raw_gb'].sum():.1f}" if not surv.empty else '-'],
        ['Total Gb (clean)', f"{surv['clean_gb'].sum():.1f}" if not surv.empty else '-'],
        ['% Unclassified', f"{run_stats['pct_unclassified']}%"],
        ['MAGs recuperados', str(n_mags)],
        ['Genes AMR', str(n_amr)],
        ['Riesgo CRITICO', str(n_crit)],
    ]
    cover_table = Table(cover_data, colWidths=[6*cm, 4*cm])
    cover_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (0,-1), C_BLUE),
        ('TEXTCOLOR', (0,0), (0,-1), C_WHITE),
        ('FONTNAME', (0,0), (-1,-1), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,-1), 10),
        ('ALIGN', (1,0), (1,-1), 'CENTER'),
        ('GRID', (0,0), (-1,-1), 0.5, C_GREY),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ]))
    story.append(cover_table)
    story.append(PageBreak())

    # ── A. ESTADISTICAS CARRERA ────────────────────────────────
    story.append(Paragraph('A. Estadisticas de la Carrera de Secuenciacion', st['SectionH']))
    story.append(Paragraph(
        f'Total reads secuenciados: <b>{run_stats["total_reads"]:,}</b>. '
        f'Reads con barcode: <b>{run_stats["classified_reads"]:,}</b>. '
        f'Reads sin barcode (unclassified): <b>{run_stats["unclassified_reads"]:,}</b> '
        f'(<b>{run_stats["pct_unclassified"]}%</b>). '
        f'{"Valor aceptable." if run_stats["pct_unclassified"]<15 else "ALERTA: >=15% de lecturas sin barcode."}'
        , st['Body']))

    # ── B. RETENCION ───────────────────────────────────────────
    story.append(Paragraph('B. Retencion de Lecturas', st['SectionH']))
    if not surv.empty:
        headers = ['Muestra', 'Raw Reads', 'Raw Gb', 'Clean Reads', 'Clean Gb', '% Ret.']
        rows_data = []
        for _, r in surv.iterrows():
            if 'unclassified' in str(r['Sample']).lower(): continue
            rows_data.append([r['Sample'], f"{r['raw_reads']:,.0f}", f"{r['raw_gb']:.2f}",
                            f"{r['clean_reads']:,.0f}", f"{r['clean_gb']:.2f}", f"{r['pct']:.1f}%"])
        story.append(make_table(headers, rows_data,
                     col_widths=[3.5*cm, 2.5*cm, 1.8*cm, 2.5*cm, 1.8*cm, 1.5*cm]))
        story.append(Paragraph(
            f'Retencion media: <b>{surv["pct"].mean():.1f}%</b> (std: {surv["pct"].std():.1f}%). '
            f'Total Gb clean: <b>{surv["clean_gb"].sum():.2f}</b>.', st['Small']))

    # ── C. CATALOGO MAGs ───────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('C. Catalogo de MAGs', st['SectionH']))
    if not checkm2.empty:
        n_hq = len(checkm2[checkm2['Q']=='HQ'])
        n_mq = len(checkm2[checkm2['Q']=='MQ'])
        n_lq = len(checkm2[checkm2['Q']=='LQ'])
        story.append(Paragraph(
            f'Total: <b>{len(checkm2)} MAGs</b> recuperados. '
            f'<font color="#28A745"><b>{n_hq} HQ</b></font>, '
            f'<font color="#E05C2A"><b>{n_mq} MQ</b></font>, '
            f'<font color="#6C757D"><b>{n_lq} LQ</b></font> (estandares MIMAG).', st['Body']))

        df = checkm2.copy()
        if not tax.empty:
            df = df.merge(tax[['Sample','Bin','Org']], on=['Sample','Bin'], how='left').fillna('')
        else:
            df['Org'] = ''

        headers = ['Muestra', 'Bin', 'Organismo', 'Comp.%', 'Cont.%', 'Q']
        rows_data = [[r['Sample'], r['Bin'][:25], r.get('Org','')[:30],
                      f"{r['Comp']:.1f}", f"{r['Cont']:.1f}", r['Q']]
                     for _, r in df.iterrows()]
        story.append(make_table(headers, rows_data,
                     col_widths=[2.8*cm, 3*cm, 4*cm, 1.5*cm, 1.5*cm, 1*cm]))

    # ── D. AMR DETALLADO ───────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('D. Resistencia Antimicrobiana en MAGs', st['SectionH']))
    if not integ.empty:
        n_plas = len(integ[integ['Location']=='PLASMID'])
        story.append(Paragraph(
            f'Total: <b>{len(integ)} genes AMR</b> detectados en <b>{integ["MAG"].nunique()} MAGs</b>. '
            f'<b>{n_plas}</b> en plasmido ({n_plas/len(integ)*100:.0f}% transferibles). '
            f'<font color="#E74C3C"><b>{n_crit} CRITICOS</b></font>.', st['Body']))

        if n_crit > 0:
            story.append(Paragraph('ALERTAS CRITICAS:', st['Alert']))
            for _, r in integ[integ['Risk_level']=='CRITICO'].iterrows():
                story.append(Paragraph(
                    f'&bull; <b>{r.get("Gene","")}</b> ({r.get("AMR_class","")}) en '
                    f'<i>{r.get("Organism","")}</i> — plasmido score {r.get("Plasmid_score",0):.2f} '
                    f'[{r.get("Sample","")}]', st['Body']))

        headers = ['Muestra', 'Organismo', 'Gen', 'Clase', 'Ubic.', 'Id%', 'Cov%', 'Riesgo']
        rows_data = []
        for _, r in integ.iterrows():
            loc = 'PLASM' if r.get('Location')=='PLASMID' else 'CROM'
            rows_data.append([r.get('Sample',''), r.get('Organism','')[:25],
                            r.get('Gene',''), r.get('AMR_class','')[:15], loc,
                            f"{r.get('Identity_pct',0):.0f}", f"{r.get('Coverage_pct',0):.0f}",
                            r.get('Risk_level','')])
        story.append(make_table(headers, rows_data,
                     col_widths=[2.5*cm, 3*cm, 1.8*cm, 2*cm, 1.2*cm, 1*cm, 1*cm, 1.3*cm]))
        story.append(Paragraph('Ver Excel adjunto para campos completos: accession, closest reference, metodo, HMM.', st['Small']))

    # ── E. CONCORDANCIA KMA-MAG ────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('E. Concordancia KMA Reads vs MAG Contigs', st['SectionH']))
    story.append(Paragraph(
        'Genes AMR detectados en reads (KMA/ResFinder) vs contigs ensamblados (AMRFinderPlus). '
        'Deteccion en ambos = confianza maxima. Solo en reads = posible reservorio no ensamblado.', st['Body']))

    if not kma.empty or not integ.empty:
        all_samples = set()
        if not kma.empty: all_samples.update(kma['Sample'].unique())
        if not integ.empty: all_samples.update(integ['Sample'].unique())
        conc_rows = []
        for s in sorted(all_samples):
            kma_genes = set(kma[kma['Sample']==s]['Gene'].unique()) if not kma.empty else set()
            mag_genes = set(integ[integ['Sample']==s]['Gene'].unique()) if not integ.empty else set()
            both = kma_genes & mag_genes
            only_reads = kma_genes - mag_genes
            only_contigs = mag_genes - kma_genes
            for g in sorted(both): conc_rows.append([s, g, 'AMBOS', 'MAXIMA'])
            for g in sorted(only_reads): conc_rows.append([s, g, 'Solo Reads', 'Moderada'])
            for g in sorted(only_contigs): conc_rows.append([s, g, 'Solo Contigs', 'Moderada'])
        if conc_rows:
            story.append(make_table(['Muestra', 'Gen', 'Deteccion', 'Confianza'], conc_rows,
                         col_widths=[3*cm, 3.5*cm, 3*cm, 2.5*cm]))

    # ── F. VIRULENCIA ──────────────────────────────────────────
    story.append(Paragraph('F. Factores de Virulencia (VFDB)', st['SectionH']))
    if not vfdb.empty:
        story.append(Paragraph(f'<b>{len(vfdb)} factores de virulencia</b> detectados por ABRicate/VFDB.', st['Body']))
        headers = ['Muestra', 'Bin', 'Gen', 'Producto', 'Id%', 'Cov%']
        rows_data = [[r['Sample'], r['Bin'][:20], r['Gene'], r['Product'][:30],
                      f"{r['Identity']:.0f}", f"{r['Coverage']:.0f}"] for _, r in vfdb.iterrows()]
        story.append(make_table(headers, rows_data,
                     col_widths=[2.5*cm, 2.5*cm, 2*cm, 4*cm, 1.2*cm, 1.2*cm]))
    else:
        story.append(Paragraph('Sin factores de virulencia detectados por VFDB.', st['Body']))

    # ── G. PLASMIDOS ───────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('G. Plasmidos y Elementos Moviles', st['SectionH']))
    if not data['plasmid_counts'].empty:
        story.append(Paragraph('Contigs clasificados como plasmidicos por geNomad:', st['Body']))
        headers = ['Muestra', 'N Plasmidos detectados']
        rows_data = [[r['Sample'], str(r['N_plasmids'])]
                     for _, r in data['plasmid_counts'].iterrows()]
        story.append(make_table(headers, rows_data, col_widths=[5*cm, 4*cm]))
    if not integ.empty:
        plas_amr = integ[integ['Location']=='PLASMID']
        if not plas_amr.empty:
            story.append(Spacer(1, 0.5*cm))
            story.append(Paragraph(f'<b>{len(plas_amr)} genes AMR</b> localizados en contigs plasmidicos:', st['Body']))
            headers = ['Muestra', 'Organismo', 'Gen', 'Clase', 'Plasmid Score', 'Riesgo']
            rows_data = [[r.get('Sample',''), r.get('Organism','')[:25], r.get('Gene',''),
                         r.get('AMR_class','')[:15], f"{r.get('Plasmid_score',0):.2f}",
                         r.get('Risk_level','')] for _, r in plas_amr.iterrows()]
            story.append(make_table(headers, rows_data,
                         col_widths=[2.5*cm, 3*cm, 2*cm, 2*cm, 2*cm, 1.5*cm]))
    story.append(Paragraph('Ver hoja "Plasmidos_MOBsuite" del Excel para tipificacion completa (replicones, movilidad, host range).', st['Small']))

    # ── H. DESGLOSE POR MUESTRA ───────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('H. Desglose por Muestra', st['SectionH']))

    samples = sorted(surv['Sample'].unique()) if not surv.empty else []
    for sample in samples:
        if 'unclassified' in sample.lower(): continue

        story.append(Paragraph(f'Muestra: {sample}', st['SubH']))

        # QC
        s_surv = surv[surv['Sample']==sample]
        if not s_surv.empty:
            r = s_surv.iloc[0]
            story.append(Paragraph(
                f'Raw: {r["raw_reads"]:,.0f} reads ({r["raw_gb"]:.2f} Gb) | '
                f'Clean: {r["clean_reads"]:,.0f} reads ({r["clean_gb"]:.2f} Gb) | '
                f'Retencion: {r["pct"]:.1f}%', st['Body']))

        # MAGs
        s_mags = checkm2[checkm2['Sample']==sample] if not checkm2.empty else pd.DataFrame()
        if not s_mags.empty:
            if not tax.empty:
                s_mags = s_mags.merge(tax[['Sample','Bin','Org']], on=['Sample','Bin'], how='left').fillna('')
            story.append(Paragraph(f'MAGs: {len(s_mags)} ({len(s_mags[s_mags["Q"]=="HQ"])} HQ, {len(s_mags[s_mags["Q"]=="MQ"])} MQ)', st['Body']))
            mag_rows = [[r['Bin'][:30], r.get('Org','')[:30], f"{r['Comp']:.0f}%", f"{r['Cont']:.0f}%", r['Q']]
                       for _, r in s_mags.iterrows()]
            story.append(make_table(['Bin', 'Organismo', 'Comp', 'Cont', 'Q'], mag_rows,
                         col_widths=[3.5*cm, 4*cm, 1.5*cm, 1.5*cm, 1*cm]))
        else:
            story.append(Paragraph('Sin MAGs recuperados.', st['Small']))

        # AMR
        s_amr = integ[integ['Sample']==sample] if not integ.empty else pd.DataFrame()
        if not s_amr.empty:
            n_p = len(s_amr[s_amr['Location']=='PLASMID'])
            story.append(Paragraph(f'AMR: {len(s_amr)} genes ({n_p} en plasmido)', st['Body']))
            amr_rows = [[r.get('Organism','')[:25], r.get('Gene',''), r.get('AMR_class','')[:15],
                        'PLASM' if r.get('Location')=='PLASMID' else 'CROM',
                        r.get('Risk_level','')] for _, r in s_amr.iterrows()]
            story.append(make_table(['Organismo', 'Gen', 'Clase', 'Ubic.', 'Riesgo'], amr_rows,
                         col_widths=[3.5*cm, 2*cm, 2.5*cm, 1.5*cm, 1.5*cm]))
        else:
            story.append(Paragraph('Sin genes AMR detectados en MAGs.', st['Small']))

        # Virulence
        s_vf = vfdb[vfdb['Sample']==sample] if not vfdb.empty else pd.DataFrame()
        if not s_vf.empty:
            story.append(Paragraph(f'Virulencia: {len(s_vf)} factores detectados', st['Body']))
            vf_rows = [[r['Gene'], r['Product'][:35], f"{r['Identity']:.0f}%"]
                      for _, r in s_vf.iterrows()]
            story.append(make_table(['Gen', 'Producto', 'Identidad'], vf_rows,
                         col_widths=[3*cm, 6*cm, 2*cm]))

        story.append(Spacer(1, 0.3*cm))
        story.append(HRFlowable(width='100%', thickness=0.5, color=C_LGREY))

    # ── I. SOFTWARE ────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('I. Software y Versiones', st['SectionH']))
    sw = [
        ['NanoPlot','1.46.2','QC reads'], ['FastQC','0.12.1','QC reads'],
        ['FastQ Screen','0.16.0','Screening'], ['Porechop ABI','0.5.1','Trimming'],
        ['Chopper','0.12.0','Filtrado'], ['Kraken2','2.17.1','Taxonomia k-mers'],
        ['Bracken','3.1','Abundancias'], ['Kaiju','1.10.1','Taxonomia proteica'],
        ['Sylph','0.9.0','Perfilado ANI'], ['KMA','1.6.8','AMR reads'],
        ['MetaFlye','2.9.6','Ensamblaje'], ['Medaka','2.2.1','Polishing'],
        ['QUAST','5.3.0','QC ensamblaje'], ['MetaBAT2','2.18','Binning'],
        ['MaxBin2','2.2.7','Binning'], ['SemiBin2','2.2.1','Binning DL'],
        ['DAS Tool','1.1.7','Refinamiento'], ['CheckM2','1.1.0','QC MAGs'],
        ['GTDB-Tk','2.7.0','Taxonomia MAGs'], ['AMRFinderPlus','4.2.7','AMR MAGs'],
        ['geNomad','1.12.0','Plasmidos/virus'], ['ABRicate','1.4.0','VFDB/CARD'],
        ['MOB-suite','3.1.9','Tipificacion plasmidos'], ['IntegronFinder','2.0rc6','Integrones'],
        ['Bakta','1.12.0','Anotacion MAGs'], ['MultiQC','1.33','Reporte QC'],
    ]
    story.append(make_table(['Herramienta', 'Version', 'Funcion'], sw,
                 col_widths=[3.5*cm, 2*cm, 6*cm]))

    # Build
    doc.build(story)
    print(f'[OK] PDF: {output_path}')


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description='Genera PDF EpiTaxMAG')
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--run-name', required=True)
    parser.add_argument('--input-dir', default=None)
    parser.add_argument('--branding', default='branding/')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    print(f'[INFO] Generando PDF para {args.run_name}')
    org_info, logo_path = load_branding(args.branding)
    run_stats = load_run_stats(args.input_dir) if args.input_dir else {
        'total_reads':0, 'classified_reads':0, 'unclassified_reads':0, 'pct_unclassified':0}
    data = load_all_data(args.results_dir)

    print(f'[INFO] Datos: {len(data["survival"])} muestras, {len(data["checkm2"])} MAGs, '
          f'{len(data["integration"])} AMR, {len(data["vfdb"])} virulencia')

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    build_pdf(data, run_stats, org_info, logo_path, args.run_name, args.output)


if __name__ == '__main__':
    main()
