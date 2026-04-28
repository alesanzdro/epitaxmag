#!/usr/bin/env python3
"""
generate_excel_report.py
Genera un Excel exhaustivo multi-hoja con TODOS los datos de EpiTaxMAG.

Hojas:
  0. Resumen del Run
  1. Retencion de Lecturas
  2. Catalogo de MAGs (taxonomia + calidad + metricas)
  3. AMR Detallado (todos los campos AMRFinderPlus)
  4. AMR por Clase (matriz organismos x clases)
  5. Concordancia KMA vs MAGs
  6. Plasmidos (geNomad + integracion)
  7. Evaluacion de Riesgo (consolidado)
  8. Virulencia y Stress (AMRFinderPlus --plus)
  9. Software y Versiones

Uso:
  python3 scripts/generate_excel_report.py \
      --results-dir results/260226_EPIM232 \
      --run-name 260226_EPIM232 \
      --input-dir /path/to/raw/fastqs \
      --output results/260226_EPIM232/28_reports/260226_EPIM232_full_report.xlsx
"""

import argparse, os, glob, warnings, gzip
from datetime import datetime
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side, numbers
from openpyxl.utils import get_column_letter
from openpyxl.utils.dataframe import dataframe_to_rows

warnings.filterwarnings('ignore')

# Colores
BLUE_FILL = PatternFill(start_color='2C5F8A', end_color='2C5F8A', fill_type='solid')
LBLUE_FILL = PatternFill(start_color='D6E8F5', end_color='D6E8F5', fill_type='solid')
GREEN_FILL = PatternFill(start_color='C8E6C9', end_color='C8E6C9', fill_type='solid')
ORANGE_FILL = PatternFill(start_color='FFE0B2', end_color='FFE0B2', fill_type='solid')
RED_FILL = PatternFill(start_color='FFCDD2', end_color='FFCDD2', fill_type='solid')
GREY_FILL = PatternFill(start_color='F5F5F5', end_color='F5F5F5', fill_type='solid')
WHITE_FONT = Font(color='FFFFFF', bold=True, size=10)
HEADER_FONT = Font(bold=True, size=10)
LINK_FONT = Font(color='0563C1', underline='single', size=9)
THIN_BORDER = Border(
    left=Side(style='thin', color='CCCCCC'),
    right=Side(style='thin', color='CCCCCC'),
    top=Side(style='thin', color='CCCCCC'),
    bottom=Side(style='thin', color='CCCCCC'))

def style_header(ws, row=1):
    for cell in ws[row]:
        cell.font = WHITE_FONT
        cell.fill = BLUE_FILL
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
        cell.border = THIN_BORDER

def style_data(ws, start_row=2):
    for row in ws.iter_rows(min_row=start_row, max_row=ws.max_row):
        for cell in row:
            cell.border = THIN_BORDER
            cell.alignment = Alignment(vertical='center', wrap_text=False)
            if cell.row % 2 == 0:
                cell.fill = GREY_FILL

def auto_width(ws, min_w=8, max_w=45):
    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        lengths = [len(str(cell.value or '')) for cell in col]
        w = min(max(max(lengths) + 2, min_w), max_w)
        ws.column_dimensions[letter].width = w

def add_conditional_fill(ws, col_letter, start_row, end_row, thresholds):
    """Apply fill based on numeric thresholds: [(val, fill), ...]"""
    for row in range(start_row, end_row + 1):
        cell = ws[f'{col_letter}{row}']
        try:
            v = float(cell.value)
            for threshold, fill in sorted(thresholds, reverse=True):
                if v >= threshold:
                    cell.fill = fill
                    break
        except (ValueError, TypeError):
            pass


# ═══════════════════════════════════════════════════════════
# DATA LOADERS (same logic as generate_mag_report.py)
# ═══════════════════════════════════════════════════════════

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
    stats = {'total_reads': 0, 'classified_reads': 0, 'unclassified_reads': 0, 'pct_unclassified': 0}
    if not input_dir or not os.path.isdir(input_dir): return stats
    for f in sorted(glob.glob(f'{input_dir}/*.fastq.gz')):
        name = os.path.basename(f).replace('.fastq.gz', '')
        try:
            n = sum(1 for line in gzip.open(f, 'rt') if line.startswith('@'))
            stats['total_reads'] += n
            if 'unclassified' in name.lower():
                stats['unclassified_reads'] = n
            else:
                stats['classified_reads'] += n
        except: pass
    if stats['total_reads'] > 0:
        stats['pct_unclassified'] = round(stats['unclassified_reads'] / stats['total_reads'] * 100, 1)
    return stats

def load_checkm2(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/23_binqc_checkm2/*/quality_report.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                comp = float(r['Completeness']); cont = float(r['Contamination'])
                qc = 'HQ' if comp >= 90 and cont <= 5 else ('MQ' if comp >= 50 and cont <= 10 else 'LQ')
                rows.append({'Sample': sample, 'Bin': r['Name'],
                             'Completeness': comp, 'Contamination': cont, 'Quality': qc})
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_taxonomy(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/24_taxonomy_gtdbtk/*/*_taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                cls = str(r.iloc[1]) if len(r) > 1 else ''
                levels = {}
                for part in cls.split(';'):
                    p = part.strip()
                    for px, lv in [('d__','Domain'),('p__','Phylum'),('c__','Class'),
                                   ('o__','Order'),('f__','Family'),('g__','Genus'),('s__','Species')]:
                        if p.startswith(px): levels[lv] = p[len(px):] or ''
                rows.append({'Sample': sample, 'Bin': str(r.iloc[0]).strip(),
                             'Full_taxonomy': cls, **levels})
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_amr_full(rd):
    """Load ALL columns from AMRFinderPlus per_bin files."""
    rows = []
    for per_bin_dir in sorted(glob.glob(f'{rd}/25_amr_mags/*/per_bin')):
        sample = os.path.basename(os.path.dirname(per_bin_dir))
        for f in sorted(glob.glob(f'{per_bin_dir}/*_amr.tsv')):
            try:
                df = pd.read_csv(f, sep='\t')
                for _, r in df.iterrows():
                    rows.append({
                        'Sample': sample,
                        'Bin': str(r.get('Name', '')).strip(),
                        'Contig': str(r.get('Contig id', '')).strip(),
                        'Start': int(r.get('Start', 0)),
                        'Stop': int(r.get('Stop', 0)),
                        'Strand': str(r.get('Strand', '')).strip(),
                        'Gene_symbol': str(r.get('Element symbol', '')).strip(),
                        'Gene_name': str(r.get('Element name', '')).strip(),
                        'Scope': str(r.get('Scope', '')).strip(),
                        'Type': str(r.get('Type', '')).strip(),
                        'Subtype': str(r.get('Subtype', '')).strip(),
                        'Class': str(r.get('Class', '')).strip(),
                        'Subclass': str(r.get('Subclass', '')).strip(),
                        'Method': str(r.get('Method', '')).strip(),
                        'Target_length': int(r.get('Target length', 0)),
                        'Ref_length': int(r.get('Reference sequence length', 0)),
                        'Pct_coverage': float(r.get('% Coverage of reference', 0)),
                        'Pct_identity': float(r.get('% Identity to reference', 0)),
                        'Alignment_length': int(r.get('Alignment length', 0)),
                        'Accession': str(r.get('Closest reference accession', '')).strip(),
                        'Closest_reference': str(r.get('Closest reference name', '')).strip(),
                        'HMM_accession': str(r.get('HMM accession', '')).strip(),
                        'HMM_description': str(r.get('HMM description', '')).strip(),
                    })
            except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_kma(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/11_amr_kma/*.res')):
        sample = os.path.basename(f).replace('.res', '')
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                template = str(r.iloc[0]).strip()
                if not template or template.startswith('#'): continue
                rows.append({
                    'Sample': sample, 'Template': template,
                    'Gene_root': template.split('_')[0] if '_' in template else template,
                    'Score': float(r.get('Score', 0)),
                    'Expected': float(r.get('Expected', 0)),
                    'Template_length': int(r.get('Template_length', 0)),
                    'Template_Identity': float(r.get('Template_Identity', 0)),
                    'Template_Coverage': float(r.get('Template_Coverage', 0)),
                    'Query_Identity': float(r.get('Query_Identity', 0)),
                    'Query_Coverage': float(r.get('Query_Coverage', 0)),
                    'Depth': float(r.get('Depth', 0)),
                    'q_value': float(r.get('q_value', 0)),
                    'p_value': float(r.get('p_value', 0)),
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_plasmids(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/26_genomad/*/*_plasmid_summary.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                row = {'Sample': sample}
                for c in df.columns:
                    row[c] = r[c]
                row['Contig_clean'] = str(r.iloc[0]).strip().split('|')[0]
                rows.append(row)
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_abricate_vfdb(rd):
    """Load ABRicate VFDB virulence results from Phase A."""
    rows = []
    for f in sorted(glob.glob(f'{rd}/28_abricate/*/*.vfdb.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            if df.empty: continue
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample,
                    'Bin': os.path.basename(f).replace('.vfdb.tsv', ''),
                    'Contig': str(r.get('SEQUENCE', '')),
                    'Gene': str(r.get('GENE', '')),
                    'Product': str(r.get('PRODUCT', '')),
                    'Pct_coverage': float(r.get('%COVERAGE', 0)),
                    'Pct_identity': float(r.get('%IDENTITY', 0)),
                    'Accession': str(r.get('ACCESSION', '')),
                    'Database': str(r.get('DATABASE', 'vfdb')),
                    'Resistance': str(r.get('RESISTANCE', '')),
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_abricate_card(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/28_abricate/*/*.card.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            if df.empty: continue
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample,
                    'Bin': os.path.basename(f).replace('.card.tsv', ''),
                    'Contig': str(r.get('SEQUENCE', '')),
                    'Gene': str(r.get('GENE', '')),
                    'Product': str(r.get('PRODUCT', '')),
                    'Pct_coverage': float(r.get('%COVERAGE', 0)),
                    'Pct_identity': float(r.get('%IDENTITY', 0)),
                    'Accession': str(r.get('ACCESSION', '')),
                    'Database': 'card',
                    'Resistance': str(r.get('RESISTANCE', '')),
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_mobsuite(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/29_mobsuite/*/mobtyper_results.txt')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample,
                    'Contig': str(r.get('sample_id', r.get('file_id', ''))),
                    'Replicon_type': str(r.get('rep_type(s)', r.get('rep_type', ''))),
                    'Relaxase_type': str(r.get('relaxase_type(s)', r.get('relaxase_type', ''))),
                    'MPF_type': str(r.get('mpf_type', '')),
                    'Orit_type': str(r.get('orit_type(s)', r.get('orit_type', ''))),
                    'Mobility': str(r.get('predicted_mobility', '')),
                    'Host_range': str(r.get('predicted_host_range_overall_rank', '')),
                    'Nearest_neighbor': str(r.get('mash_nearest_neighbor', '')),
                    'Mash_distance': r.get('mash_neighbor_distance', ''),
                    'GC': r.get('gc', ''),
                    'Size': r.get('size', ''),
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_integronfinder(rd):
    rows = []
    for integ_dir in sorted(glob.glob(f'{rd}/30_integronfinder/*')):
        sample = os.path.basename(integ_dir)
        for f in glob.glob(f'{integ_dir}/**/*.integrons', recursive=True):
            try:
                df = pd.read_csv(f, sep='\t', comment='#')
                if df.empty: continue
                for _, r in df.iterrows():
                    rows.append({'Sample': sample, **{c: r[c] for c in df.columns}})
            except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_bakta_summary(rd):
    rows = []
    for f in sorted(glob.glob(f'{rd}/31_bakta/*/*/*.tsv')):
        bin_name = os.path.basename(os.path.dirname(f))
        sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            if df.empty: continue
            n_cds = len(df[df.iloc[:, 1] == 'cds']) if len(df.columns) > 1 else 0
            n_rrna = len(df[df.iloc[:, 1] == 'rRNA']) if len(df.columns) > 1 else 0
            n_trna = len(df[df.iloc[:, 1] == 'tRNA']) if len(df.columns) > 1 else 0
            n_crispr = len(df[df.iloc[:, 1] == 'crispr']) if len(df.columns) > 1 else 0
            rows.append({
                'Sample': sample, 'Bin': bin_name,
                'Total_features': len(df), 'CDS': n_cds,
                'rRNA': n_rrna, 'tRNA': n_trna, 'CRISPR': n_crispr
            })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def load_integration(rd):
    dfs = []
    for f in sorted(glob.glob(f'{rd}/27_amr_pathogen_integration/*_amr_pathogen.tsv')):
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty: dfs.append(df)
        except: pass
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


# ═══════════════════════════════════════════════════════════
# EXCEL GENERATION
# ═══════════════════════════════════════════════════════════

def write_df_to_sheet(ws, df, start_row=1):
    """Write DataFrame to worksheet with proper formatting."""
    # Headers
    for col_idx, col_name in enumerate(df.columns, 1):
        ws.cell(row=start_row, column=col_idx, value=col_name)
    style_header(ws, start_row)

    # Data
    for row_idx, (_, row) in enumerate(df.iterrows(), start_row + 1):
        for col_idx, val in enumerate(row, 1):
            cell = ws.cell(row=row_idx, column=col_idx)
            if pd.isna(val):
                cell.value = ''
            elif isinstance(val, float):
                cell.value = round(val, 2)
                cell.number_format = '0.00'
            else:
                cell.value = val
            cell.border = THIN_BORDER
            if row_idx % 2 == 0:
                cell.fill = GREY_FILL

    auto_width(ws)


def create_workbook(data, run_name, run_stats):
    wb = Workbook()
    date_str = datetime.now().strftime('%Y-%m-%d %H:%M')

    # ── HOJA 0: Resumen ──────────────────────────────────────
    ws = wb.active
    ws.title = 'Resumen'
    ws.sheet_properties.tabColor = '2C5F8A'

    info = [
        ('EpiTaxMAG Report', ''),
        ('Run', run_name),
        ('Fecha', date_str),
        ('', ''),
        ('ESTADISTICAS DE CARRERA', ''),
        ('Total reads secuenciados', f"{run_stats['total_reads']:,}"),
        ('Reads clasificados (barcode)', f"{run_stats['classified_reads']:,}"),
        ('Reads unclassified', f"{run_stats['unclassified_reads']:,}"),
        ('% Unclassified', f"{run_stats['pct_unclassified']}%"),
        ('', ''),
        ('DATOS PROCESADOS', ''),
    ]
    surv = data['survival']
    if not surv.empty:
        info.extend([
            ('Muestras', len(surv)),
            ('Gb raw total', f"{surv['Raw_Gb'].sum():.2f}"),
            ('Gb clean total', f"{surv['Clean_Gb'].sum():.2f}"),
            ('Retencion media', f"{surv['Retention_%'].mean():.1f}%"),
            ('Retencion std', f"{surv['Retention_%'].std():.1f}%"),
        ])
    checkm2 = data['checkm2']
    if not checkm2.empty:
        info.extend([
            ('', ''),
            ('MAGs', ''),
            ('Total MAGs', len(checkm2)),
            ('High Quality (HQ)', len(checkm2[checkm2['Quality'] == 'HQ'])),
            ('Medium Quality (MQ)', len(checkm2[checkm2['Quality'] == 'MQ'])),
            ('Low Quality (LQ)', len(checkm2[checkm2['Quality'] == 'LQ'])),
        ])
    integ = data['integration']
    if not integ.empty:
        info.extend([
            ('', ''),
            ('AMR', ''),
            ('Total genes AMR', len(integ)),
            ('En plasmido', len(integ[integ['Location'] == 'PLASMID'])),
            ('Riesgo CRITICO', len(integ[integ['Risk_level'] == 'CRITICO'])),
            ('Riesgo ALTO', len(integ[integ['Risk_level'] == 'ALTO'])),
        ])

    for row_idx, (k, v) in enumerate(info, 1):
        ws.cell(row=row_idx, column=1, value=k).font = Font(bold=True, size=11 if row_idx <= 3 else 10)
        ws.cell(row=row_idx, column=2, value=v)
    ws.column_dimensions['A'].width = 30
    ws.column_dimensions['B'].width = 25

    # ── HOJA 1: Retencion ────────────────────────────────────
    if not surv.empty:
        ws1 = wb.create_sheet('Retencion')
        ws1.sheet_properties.tabColor = '28A745'
        write_df_to_sheet(ws1, surv)

    # ── HOJA 2: Catalogo MAGs ────────────────────────────────
    tax = data['taxonomy']
    if not checkm2.empty:
        ws2 = wb.create_sheet('MAGs')
        ws2.sheet_properties.tabColor = '5B9ABF'
        if not tax.empty:
            mags = checkm2.merge(tax, on=['Sample', 'Bin'], how='left')
        else:
            mags = checkm2.copy()
        mags = mags.fillna('')
        write_df_to_sheet(ws2, mags)

    # ── HOJA 3: AMR Detallado (TODOS los campos) ────────────
    amr_full = data['amr_full']
    if not amr_full.empty:
        ws3 = wb.create_sheet('AMR Detallado')
        ws3.sheet_properties.tabColor = 'E74C3C'
        # Add taxonomy info
        if not tax.empty:
            amr_enriched = amr_full.merge(
                tax[['Sample', 'Bin', 'Genus', 'Species', 'Full_taxonomy']].drop_duplicates(),
                on=['Sample', 'Bin'], how='left')
        else:
            amr_enriched = amr_full.copy()
        # Add plasmid info
        plasmids = data['plasmids']
        if not plasmids.empty:
            plas_set = set(plasmids['Contig_clean'].values)
            plas_scores = dict(zip(plasmids['Contig_clean'],
                                   plasmids.get('plasmid_score', plasmids.iloc[:, 1])))
            amr_enriched['Location'] = amr_enriched['Contig'].apply(
                lambda c: 'PLASMID' if c in plas_set else 'CHROMOSOME')
            amr_enriched['Plasmid_score'] = amr_enriched['Contig'].apply(
                lambda c: round(plas_scores.get(c, 0), 3))
        else:
            amr_enriched['Location'] = 'UNKNOWN'
            amr_enriched['Plasmid_score'] = 0

        # Add NCBI link column
        amr_enriched['NCBI_link'] = amr_enriched['Accession'].apply(
            lambda a: f'https://www.ncbi.nlm.nih.gov/protein/{a}' if a and a != 'NA' else '')
        amr_enriched['CARD_link'] = amr_enriched['Gene_symbol'].apply(
            lambda g: f'https://card.mcmaster.ca/ontology/query?query={g}' if g else '')

        amr_enriched = amr_enriched.fillna('')
        write_df_to_sheet(ws3, amr_enriched)

        # Add hyperlinks
        acc_col = list(amr_enriched.columns).index('NCBI_link') + 1
        for row_idx in range(2, len(amr_enriched) + 2):
            cell = ws3.cell(row=row_idx, column=acc_col)
            if cell.value and cell.value.startswith('http'):
                cell.hyperlink = cell.value
                cell.font = LINK_FONT

    # ── HOJA 4: AMR por Clase (matriz) ──────────────────────
    if not integ.empty:
        ws4 = wb.create_sheet('AMR x Clase')
        ws4.sheet_properties.tabColor = 'E05C2A'
        pivot = integ.groupby(['Organism', 'AMR_class']).size().unstack(fill_value=0)
        pivot = pivot.reset_index()
        write_df_to_sheet(ws4, pivot)

    # ── HOJA 5: KMA Reads (TODOS los campos) ────────────────
    kma = data['kma']
    if not kma.empty:
        ws5 = wb.create_sheet('KMA Reads')
        ws5.sheet_properties.tabColor = 'F39C12'
        write_df_to_sheet(ws5, kma)

    # ── HOJA 6: Concordancia ─────────────────────────────────
    if not integ.empty and not kma.empty:
        ws6 = wb.create_sheet('Concordancia')
        ws6.sheet_properties.tabColor = '1ABC9C'
        # Build concordance
        conc_rows = []
        for sample in sorted(set(kma['Sample'].unique()) | set(amr_full['Sample'].unique() if not amr_full.empty else [])):
            kma_genes = set(kma[kma['Sample'] == sample]['Gene_root'].unique()) if not kma.empty else set()
            mag_genes = set(amr_full[amr_full['Sample'] == sample]['Gene_symbol'].unique()) if not amr_full.empty else set()
            for g in sorted(kma_genes | mag_genes):
                in_kma = g in kma_genes
                in_mag = g in mag_genes
                source = 'AMBOS' if in_kma and in_mag else ('SOLO_READS' if in_kma else 'SOLO_CONTIGS')
                kma_row = kma[(kma['Sample'] == sample) & (kma['Gene_root'] == g)]
                depth = kma_row['Depth'].max() if not kma_row.empty else 0
                conc_rows.append({'Sample': sample, 'Gene': g, 'Deteccion': source,
                                  'KMA_Depth': round(depth, 1),
                                  'Confianza': 'MAXIMA' if source == 'AMBOS' else 'MODERADA'})
        conc_df = pd.DataFrame(conc_rows)
        write_df_to_sheet(ws6, conc_df)

    # ── HOJA 7: Plasmidos ────────────────────────────────────
    plasmids = data['plasmids']
    if not plasmids.empty:
        ws7 = wb.create_sheet('Plasmidos')
        ws7.sheet_properties.tabColor = '9B59B6'
        write_df_to_sheet(ws7, plasmids)

    # ── HOJA 8: Riesgo Consolidado ───────────────────────────
    if not integ.empty:
        ws8 = wb.create_sheet('Riesgo')
        ws8.sheet_properties.tabColor = 'E74C3C'
        risk = integ.sort_values('Risk_level', key=lambda x: x.map(
            {'CRITICO': 0, 'ALTO': 1, 'MEDIO': 2, 'BAJO': 3}))
        write_df_to_sheet(ws8, risk)

    # ── HOJA 9: Virulencia y Stress ──────────────────────────
    if not amr_full.empty:
        vir_stress = amr_full[amr_full['Type'].isin(['VIRULENCE', 'STRESS'])].copy()
        if not vir_stress.empty:
            ws9 = wb.create_sheet('Virulencia_Stress')
            ws9.sheet_properties.tabColor = 'FF5722'
            write_df_to_sheet(ws9, vir_stress)

    # ── HOJAS FASE A (si existen) ───────────────────────────
    vfdb = data.get('vfdb', pd.DataFrame())
    if not vfdb.empty:
        ws_vf = wb.create_sheet('Virulencia_VFDB')
        ws_vf.sheet_properties.tabColor = 'FF5722'
        write_df_to_sheet(ws_vf, vfdb)

    card = data.get('card', pd.DataFrame())
    if not card.empty:
        ws_cd = wb.create_sheet('AMR_CARD')
        ws_cd.sheet_properties.tabColor = 'D32F2F'
        write_df_to_sheet(ws_cd, card)

    mob = data.get('mobsuite', pd.DataFrame())
    if not mob.empty:
        ws_mob = wb.create_sheet('Plasmidos_MOBsuite')
        ws_mob.sheet_properties.tabColor = '7B1FA2'
        write_df_to_sheet(ws_mob, mob)

    integ_finder = data.get('integronfinder', pd.DataFrame())
    if not integ_finder.empty:
        ws_if = wb.create_sheet('Integrones')
        ws_if.sheet_properties.tabColor = '00695C'
        write_df_to_sheet(ws_if, integ_finder)

    bakta = data.get('bakta', pd.DataFrame())
    if not bakta.empty:
        ws_bk = wb.create_sheet('Bakta_Anotacion')
        ws_bk.sheet_properties.tabColor = '0277BD'
        write_df_to_sheet(ws_bk, bakta)

    # ── HOJA FINAL: Software ────────────────────────────────
    ws10 = wb.create_sheet('Software')
    ws10.sheet_properties.tabColor = '607D8B'
    sw = [
        ('Herramienta', 'Version', 'Fase', 'Funcion'),
        ('NanoPlot', '1.46.2', 'TAX', 'QC reads'),
        ('FastQC', '0.12.1', 'TAX', 'QC reads'),
        ('FastQ Screen', '0.16.0', 'TAX', 'Screening contaminacion'),
        ('Porechop ABI', '0.5.1', 'TAX', 'Trimming adaptadores'),
        ('Chopper', '0.12.0', 'TAX', 'Filtrado Q + longitud'),
        ('minimap2', '2.30', 'TAX/MAG', 'Mapeo'),
        ('samtools', '1.21', 'TAX/MAG', 'BAM processing'),
        ('Kraken2', '2.17.1', 'TAX', 'Taxonomia k-mers'),
        ('Bracken', '3.1', 'TAX', 'Abundancias'),
        ('Kaiju', '1.10.1', 'TAX', 'Taxonomia proteica'),
        ('Sylph', '0.9.0', 'TAX', 'Perfilado ANI'),
        ('KMA', '1.6.8', 'TAX', 'AMR desde reads (ResFinder)'),
        ('MultiQC', '1.33', 'TAX', 'Reporte integrado'),
        ('MetaFlye', '2.9.6', 'MAG', 'Ensamblaje metagenomico'),
        ('Medaka', '2.2.1', 'MAG', 'Polishing'),
        ('QUAST', '5.3.0', 'MAG', 'QC ensamblaje'),
        ('MetaBAT2', '2.18', 'MAG', 'Binning'),
        ('MaxBin2', '2.2.7', 'MAG', 'Binning'),
        ('SemiBin2', '2.2.1', 'MAG', 'Binning deep learning'),
        ('DAS Tool', '1.1.7', 'MAG', 'Refinamiento bins'),
        ('CheckM2', '1.1.0', 'MAG', 'QC MAGs'),
        ('GTDB-Tk', '2.7.0', 'MAG', 'Taxonomia MAGs (r232)'),
        ('AMRFinderPlus', '4.2.7', 'MAG', 'AMR + virulencia + stress'),
        ('geNomad', '1.12.0', 'MAG', 'Plasmidos y virus'),
        ('ABRicate', '1.4.0', 'Fase A', 'Virulencia (VFDB) + CARD + PlasmidFinder'),
        ('MOB-suite', '3.1.9', 'Fase A', 'Tipificacion plasmidos (replicones, movilidad)'),
        ('IntegronFinder', '2.0.6', 'Fase A', 'Deteccion de integrones'),
        ('Bakta', '1.12.0', 'Fase A', 'Anotacion funcional MAGs'),
    ]
    for row_idx, row_data in enumerate(sw, 1):
        for col_idx, val in enumerate(row_data, 1):
            ws10.cell(row=row_idx, column=col_idx, value=val)
    style_header(ws10)
    style_data(ws10)
    auto_width(ws10)

    return wb


def main():
    parser = argparse.ArgumentParser(description='Genera Excel exhaustivo EpiTaxMAG')
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--run-name', required=True)
    parser.add_argument('--input-dir', default=None)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    rd = args.results_dir
    print(f'[INFO] Generando Excel para {args.run_name}')

    run_stats = load_run_stats(args.input_dir) if args.input_dir else {
        'total_reads': 0, 'classified_reads': 0, 'unclassified_reads': 0, 'pct_unclassified': 0}

    qc_raw = load_qc_csv(f'{rd}/01_qc_raw/qc_raw_metrics.csv')
    qc_filt = load_qc_csv(f'{rd}/06_qc_filtered/qc_filtered_metrics.csv')

    # Survival
    survival = pd.DataFrame()
    if not qc_raw.empty and not qc_filt.empty:
        raw = qc_raw[['Sample', 'total_reads', 'mean_length']].copy()
        raw.columns = ['Sample', 'Raw_reads', 'Raw_mean_len']
        raw['Raw_Gb'] = raw['Raw_reads'] * raw['Raw_mean_len'] / 1e9
        filt = qc_filt[['Sample', 'total_reads', 'mean_length']].copy()
        filt.columns = ['Sample', 'Clean_reads', 'Clean_mean_len']
        filt['Clean_Gb'] = filt['Clean_reads'] * filt['Clean_mean_len'] / 1e9
        survival = pd.merge(raw, filt, on='Sample', how='outer').fillna(0)
        survival['Pct_retention'] = np.where(
            survival['Raw_reads'] > 0,
            (survival['Clean_reads'] / survival['Raw_reads'] * 100).round(1), 0)
        survival = survival[['Sample', 'Raw_reads', 'Raw_Gb', 'Clean_reads', 'Clean_Gb', 'Pct_retention']]
        survival = survival.sort_values('Sample').reset_index(drop=True)
        # Rename for clarity
        survival.columns = ['Sample', 'Raw_reads', 'Raw_Gb', 'Clean_reads', 'Clean_Gb', 'Retention_%']

    data = {
        'survival': survival,
        'checkm2': load_checkm2(rd),
        'taxonomy': load_taxonomy(rd),
        'amr_full': load_amr_full(rd),
        'kma': load_kma(rd),
        'plasmids': load_plasmids(rd),
        'integration': load_integration(rd),
        # Fase A
        'vfdb': load_abricate_vfdb(rd),
        'card': load_abricate_card(rd),
        'mobsuite': load_mobsuite(rd),
        'integronfinder': load_integronfinder(rd),
        'bakta': load_bakta_summary(rd),
    }

    n_vfdb = len(data['vfdb'])
    n_mob = len(data['mobsuite'])
    print(f'[INFO] Datos: {len(survival)} muestras, {len(data["checkm2"])} MAGs, '
          f'{len(data["amr_full"])} AMR, {len(data["kma"])} KMA, {len(data["plasmids"])} plasmidos')
    print(f'[INFO] Fase A: {n_vfdb} virulencia VFDB, {n_mob} plasmidos MOB-suite, '
          f'{len(data["integronfinder"])} integrones, {len(data["bakta"])} MAGs anotados')

    wb = create_workbook(data, args.run_name, run_stats)

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    wb.save(args.output)
    print(f'[OK] Excel: {args.output}')
    print(f'[OK] Hojas: {", ".join(wb.sheetnames)}')


if __name__ == '__main__':
    main()
