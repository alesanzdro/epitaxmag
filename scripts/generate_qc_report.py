#!/usr/bin/env python3
"""
generate_qc_report.py
Genera un informe PDF de calidad Nanopore a partir del CSV de QC.
Uso: python3 generate_qc_report.py <qc_file.csv> [output.pdf]
"""

import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import io
from datetime import datetime

from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, PageBreak, Image, HRFlowable)

# ─── Paleta FISABIO ──────────────────────────────────────────────────────────
COL_BLUE    = colors.HexColor('#2C5F8A')
COL_LBLUE   = colors.HexColor('#5B9ABF')
COL_GREY    = colors.HexColor('#6C757D')
COL_LGREY   = colors.HexColor('#F0F4F7')
COL_WHITE   = colors.white
COL_GREEN   = colors.HexColor('#28A745')
COL_RED     = colors.HexColor('#DC3545')

BAR_COLOR   = '#5B9ABF'
LINE_MEAN   = '#E05C2A'
LINE_N50    = '#2C8A4A'

def parse_float(val):
    """Parsea flotantes con coma decimal (formato europeo)."""
    if isinstance(val, (int, float)):
        return float(val)
    return float(str(val).replace(',', '.'))

def load_data(csv_path):
    # Auto-detectar separador (tab o coma)
    with open(csv_path, 'r') as f:
        header = f.readline()
    sep = '\t' if '\t' in header else (';' if ';' in header else ',')
    
    df = pd.read_csv(csv_path, sep=sep)
    numeric_cols = ['mean_length', 'median_length', 'gc_percentage',
                    'dorado_mean_q', 'custom_mean_q', 'n50',
                    'total_reads', 'min_length', 'max_length']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = df[col].apply(parse_float)
    df = df.sort_values('Sample').reset_index(drop=True)
    return df

def fig_to_image(fig, dpi=150):
    """Convierte figura matplotlib a imagen para ReportLab."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    buf.seek(0)
    plt.close(fig)
    return buf

def make_piechart(df):
    """
    Donut chart central pequeño con líneas radiales hacia etiquetas externas.
    Estilo similar a pycoQC / ggrepel — sin leyenda, etiquetas directas.
    """
    df_plot = df[~df['Sample'].str.lower().str.startswith('negativo')].copy()
    n = len(df_plot)
    reads = df_plot['total_reads'].values
    labels = df_plot['Sample'].values
    total  = reads.sum()

    # Paleta garantizada para N colores (hasta 128)
    def build_palette(n):
        palette = []
        cmaps = ['tab20', 'tab20b', 'tab20c', 'Set3', 'Set2', 'Paired', 'Dark2', 'Accent']
        for cm_name in cmaps:
            cmap = matplotlib.colormaps.get_cmap(cm_name)
            nc = cmap.N if hasattr(cmap, 'N') else 20
            for i in range(nc):
                palette.append(cmap(i / nc))
        # Si aún faltan colores, rellenar con HSV distribuido
        if len(palette) < n:
            for i in range(n - len(palette)):
                hue = (i * 0.618033988749895) % 1.0  # golden ratio spacing
                import colorsys
                rgb = colorsys.hsv_to_rgb(hue, 0.65, 0.85)
                palette.append((*rgb, 1.0))
        return palette[:n]

    pie_colors = build_palette(n)

    fig, ax = plt.subplots(figsize=(14, 14), facecolor='white')
    ax.set_aspect('equal')
    ax.set_xlim(-2.8, 2.8)
    ax.set_ylim(-2.8, 2.8)
    ax.axis('off')

    # Radio del donut
    r_inner = 0.55
    r_outer = 1.0

    # Dibujar wedges manualmente para control total
    angles = np.cumsum([0] + list(reads / total * 360))
    start_angles = angles[:-1]
    end_angles   = angles[1:]

    for i in range(n):
        theta1 = start_angles[i]
        theta2 = end_angles[i]
        wedge = mpatches.Wedge(
            center=(0, 0),
            r=r_outer,
            theta1=theta1, theta2=theta2,
            width=r_outer - r_inner,
            facecolor=pie_colors[i],
            edgecolor='white',
            linewidth=0.5
        )
        ax.add_patch(wedge)

    # Texto central
    ax.text(0, 0.08, 'Total', ha='center', va='center',
            fontsize=10, fontweight='bold', color='#2C5F8A')
    ax.text(0, -0.15, f'{total:,.0f}', ha='center', va='center',
            fontsize=9, color='#444444')

    # ── Etiquetas con líneas radiales ────────────────────────────────────────
    # Radio donde empieza la línea (exterior del donut)
    r_line_start = 1.05
    # Radio donde termina la línea (zona de texto)
    r_label      = 2.55

    # Calcular ángulo medio de cada wedge y posición
    mid_angles = [(start_angles[i] + end_angles[i]) / 2 for i in range(n)]

    # Convertir a radianes y calcular coords
    label_positions = []
    for i, mid_deg in enumerate(mid_angles):
        mid_rad = np.deg2rad(mid_deg)
        # Punto de arranque (borde exterior del donut)
        x0 = r_line_start * np.cos(mid_rad)
        y0 = r_line_start * np.sin(mid_rad)
        # Punto intermedio (codo de la línea)
        r_mid = r_outer + 0.38
        xm = r_mid * np.cos(mid_rad)
        ym = r_mid * np.sin(mid_rad)
        # Punto final (etiqueta)
        xl = r_label * np.cos(mid_rad)
        yl = r_label * np.sin(mid_rad)
        label_positions.append((x0, y0, xm, ym, xl, yl, mid_deg))

    for i, (x0, y0, xm, ym, xl, yl, mid_deg) in enumerate(label_positions):
        pct = reads[i] / total * 100

        # Solo dibujar etiqueta si el segmento es visible (>0.3%)
        if pct < 0.3:
            continue

        # Línea codo: punto del wedge → codo → texto
        ha = 'left' if xl >= 0 else 'right'
        # Punto de fin de línea horizontal
        x_end = (r_label + 0.05) * np.sign(xl) if xl != 0 else 0

        ax.annotate(
            '',
            xy=(xl, yl),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle='-',
                color=pie_colors[i],
                lw=0.7,
                connectionstyle='arc3,rad=0.0'
            )
        )

        # Texto: nombre de muestra + porcentaje
        fontsize = 5.0 if n > 60 else (5.5 if n > 40 else 6.5)
        label_text = f'{labels[i]}\n{pct:.1f}%'
        ax.text(xl, yl, label_text,
                ha=ha, va='center',
                fontsize=fontsize,
                color='#222222',
                linespacing=1.2)

    ax.set_title('Distribución de lecturas por muestra',
                 fontsize=14, fontweight='bold', color='#2C5F8A',
                 y=0.97)

    return fig

def make_barline_chart(df):
    """
    Gráfico horizontal: barras (total_reads) + líneas mean_length y N50.
    Barras anchas, letra más grande, márgenes ajustados.
    """
    n = len(df)
    # Altura dinámica muy ajustada para aprovechar espacio
    fig_h = max(14, n * 0.22)
    
    fig, ax1 = plt.subplots(figsize=(15, fig_h), facecolor='white')
    # Márgenes muy ajustados
    fig.subplots_adjust(left=0.22, right=0.96, top=0.96, bottom=0.04)

    y = np.arange(n)
    bar_h = 0.75  # Barras más anchas = aspecto más moderno

    # Barras horizontales con gradiente de color por valor
    norm_reads = df['total_reads'] / df['total_reads'].max()
    bar_colors = [plt.cm.Blues(0.4 + 0.5 * v) for v in norm_reads]
    
    bars = ax1.barh(y, df['total_reads'], height=bar_h,
                    color=bar_colors, alpha=0.88,
                    edgecolor='white', linewidth=0.4,
                    label='Total reads')

    ax1.set_xlabel('Total reads', fontsize=10, color='#2C5F8A', labelpad=6)
    ax1.set_yticks(y)
    ax1.set_yticklabels(df['Sample'], fontsize=8.5, fontfamily='monospace')
    ax1.tick_params(axis='x', colors='#2C5F8A', labelsize=8.5)
    ax1.tick_params(axis='y', length=0, pad=4)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_color('#CCCCCC')
    ax1.invert_yaxis()
    ax1.set_axisbelow(True)
    ax1.grid(axis='x', linestyle='--', alpha=0.25, color='#AAAAAA', linewidth=0.7)

    # Eje secundario — Mean length y N50
    ax2 = ax1.twiny()
    ax2.plot(df['mean_length'], y, 'o-',
             color=LINE_MEAN, linewidth=1.5, markersize=4,
             label='Mean Length', alpha=0.92, zorder=5)
    ax2.plot(df['n50'], y, 's--',
             color=LINE_N50, linewidth=1.5, markersize=4,
             label='N50', alpha=0.92, zorder=5)

    ax2.set_xlabel('Longitud (bp)', fontsize=10, color='#555555', labelpad=6)
    ax2.tick_params(axis='x', colors='#555555', labelsize=8.5)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_color('#CCCCCC')

    # Leyenda combinada
    patch_reads = mpatches.Patch(color=plt.cm.Blues(0.7), alpha=0.88, label='Total reads')
    line_mean   = plt.Line2D([0],[0], color=LINE_MEAN, marker='o', markersize=5, label='Mean Length')
    line_n50    = plt.Line2D([0],[0], color=LINE_N50, marker='s', linestyle='--', markersize=5, label='N50')
    ax1.legend(handles=[patch_reads, line_mean, line_n50],
               loc='lower right', fontsize=9, framealpha=0.9,
               edgecolor='#CCCCCC')

    ax1.set_title('Total reads, Mean Length y N50 por muestra',
                  fontsize=13, fontweight='bold', color='#2C5F8A', pad=10)

    return fig

def make_summary_table(df):
    """Calcula estadísticas resumen."""
    cols = ['total_reads', 'mean_length', 'n50', 'gc_percentage', 'dorado_mean_q']
    labels = {
        'total_reads':    'Total reads',
        'mean_length':    'Mean length (bp)',
        'n50':            'N50 (bp)',
        'gc_percentage':  'GC (%)',
        'dorado_mean_q':  'Dorado mean Q'
    }
    rows = []
    for col in cols:
        if col in df.columns:
            vals = df[col]
            rows.append([
                labels[col],
                f'{vals.mean():,.2f}',
                f'{vals.median():,.2f}',
                f'{vals.std():,.2f}',
                f'{vals.min():,.2f}',
                f'{vals.max():,.2f}'
            ])
    return rows

def build_pdf(csv_path, output_path):
    run_name = os.path.splitext(os.path.basename(csv_path))[0].replace('qc_', '')
    df = load_data(csv_path)
    n_samples = len(df)
    date_str = datetime.now().strftime('%d/%m/%Y %H:%M')

    doc = SimpleDocTemplate(
        output_path,
        pagesize=A4,
        rightMargin=1.5*cm, leftMargin=1.5*cm,
        topMargin=1.8*cm, bottomMargin=1.8*cm,
        title=f'QC Report - {run_name}',
        author='FISABIO Bioinformatics'
    )

    styles = getSampleStyleSheet()
    
    st_title = ParagraphStyle('Title2', parent=styles['Title'],
                               fontSize=18, textColor=COL_BLUE,
                               spaceAfter=4, alignment=TA_CENTER)
    st_sub = ParagraphStyle('Sub', parent=styles['Normal'],
                             fontSize=9, textColor=COL_GREY,
                             alignment=TA_CENTER, spaceAfter=8)
    st_h1 = ParagraphStyle('H1', parent=styles['Heading1'],
                            fontSize=12, textColor=COL_BLUE,
                            spaceBefore=10, spaceAfter=6)
    st_normal = ParagraphStyle('N', parent=styles['Normal'],
                                fontSize=8, textColor=colors.black)

    story = []

    # ── PÁGINA 1: Resumen + Tabla completa ──────────────────────────────────

    story.append(Paragraph(f'QC Report — {run_name}', st_title))
    story.append(Paragraph(
        f'Generado el {date_str} &nbsp;|&nbsp; {n_samples} muestras &nbsp;|&nbsp; '
        f'Pipeline: Dorado 1.3.1 · duplex · sup@v5.2.0',
        st_sub))
    story.append(HRFlowable(width='100%', thickness=1.5, color=COL_BLUE, spaceAfter=10))

    # Tabla resumen estadístico
    story.append(Paragraph('Estadísticas globales', st_h1))
    
    summary_data = [['Métrica', 'Media', 'Mediana', 'Desv. Est.', 'Mínimo', 'Máximo']]
    summary_data += make_summary_table(df)
    
    t_sum = Table(summary_data, colWidths=[4.5*cm, 2.8*cm, 2.8*cm, 2.8*cm, 2.8*cm, 2.8*cm])
    t_sum.setStyle(TableStyle([
        ('BACKGROUND',    (0, 0), (-1, 0),  COL_BLUE),
        ('TEXTCOLOR',     (0, 0), (-1, 0),  COL_WHITE),
        ('FONTNAME',      (0, 0), (-1, 0),  'Helvetica-Bold'),
        ('FONTSIZE',      (0, 0), (-1, 0),  9),
        ('ALIGN',         (0, 0), (-1, -1), 'CENTER'),
        ('VALIGN',        (0, 0), (-1, -1), 'MIDDLE'),
        ('ROWBACKGROUNDS',(0, 1), (-1, -1), [COL_LGREY, COL_WHITE]),
        ('FONTSIZE',      (0, 1), (-1, -1), 8),
        ('FONTNAME',      (0, 1), (0, -1),  'Helvetica-Bold'),
        ('GRID',          (0, 0), (-1, -1), 0.3, COL_GREY),
        ('TOPPADDING',    (0, 0), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('LEFTPADDING',   (0, 0), (-1, -1), 6),
        ('RIGHTPADDING',  (0, 0), (-1, -1), 6),
    ]))
    story.append(t_sum)
    story.append(Spacer(1, 14))

    # Tabla detalle por muestra
    story.append(Paragraph('Datos por muestra', st_h1))

    display_cols = ['Sample', 'total_reads', 'max_length', 'mean_length',
                    'median_length', 'n50', 'gc_percentage', 'dorado_mean_q']
    col_headers  = ['Muestra', 'Total reads', 'Max len', 'Mean len',
                    'Median len', 'N50', 'GC %', 'Q score']
    col_widths   = [4.8*cm, 1.9*cm, 1.6*cm, 1.9*cm, 1.9*cm, 1.6*cm, 1.6*cm, 1.7*cm]

    table_data = [col_headers]
    for _, row in df.iterrows():
        table_data.append([
            str(row['Sample']),
            f"{int(row['total_reads']):,}",
            f"{int(row['max_length']):,}",
            f"{row['mean_length']:,.1f}",
            f"{row['median_length']:,.1f}",
            f"{int(row['n50']):,}",
            f"{row['gc_percentage']:.2f}",
            f"{row['dorado_mean_q']:.2f}"
        ])

    t_detail = Table(table_data, colWidths=col_widths, repeatRows=1)
    
    # Colorear celdas de Q score
    q_styles = []
    for i, (_, row) in enumerate(df.iterrows(), start=1):
        q = row['dorado_mean_q']
        if q >= 20:
            q_styles.append(('BACKGROUND', (7, i), (7, i), colors.HexColor('#C8E6C9')))
        elif q >= 15:
            q_styles.append(('BACKGROUND', (7, i), (7, i), colors.HexColor('#FFF9C4')))
        else:
            q_styles.append(('BACKGROUND', (7, i), (7, i), colors.HexColor('#FFCDD2')))

    t_detail.setStyle(TableStyle([
        ('BACKGROUND',    (0, 0), (-1, 0),  COL_BLUE),
        ('TEXTCOLOR',     (0, 0), (-1, 0),  COL_WHITE),
        ('FONTNAME',      (0, 0), (-1, 0),  'Helvetica-Bold'),
        ('FONTSIZE',      (0, 0), (-1, 0),  7.5),
        ('ALIGN',         (1, 0), (-1, -1), 'RIGHT'),
        ('ALIGN',         (0, 0), (0, -1),  'LEFT'),
        ('VALIGN',        (0, 0), (-1, -1), 'MIDDLE'),
        ('ROWBACKGROUNDS',(0, 1), (-1, -1), [COL_LGREY, COL_WHITE]),
        ('FONTSIZE',      (0, 1), (-1, -1), 7),
        ('GRID',          (0, 0), (-1, -1), 0.25, colors.HexColor('#CCCCCC')),
        ('TOPPADDING',    (0, 0), (-1, -1), 2.5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 2.5),
        ('LEFTPADDING',   (0, 0), (-1, -1), 4),
        ('RIGHTPADDING',  (0, 0), (-1, -1), 4),
    ] + q_styles))
    story.append(t_detail)
    story.append(Spacer(1, 6))
    
    # Leyenda colores Q
    legend_text = (
        '<font color="#388E3C">■</font> Q≥20 (excelente) &nbsp;&nbsp;'
        '<font color="#F9A825">■</font> Q≥15 (aceptable) &nbsp;&nbsp;'
        '<font color="#C62828">■</font> Q&lt;15 (bajo)'
    )
    story.append(Paragraph(legend_text, ParagraphStyle('leg', parent=styles['Normal'],
                                                        fontSize=7, textColor=COL_GREY)))

    # ── PÁGINA 2: Piechart ───────────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('Distribución de lecturas por muestra', st_h1))
    story.append(Spacer(1, 6))

    fig_pie = make_piechart(df)
    buf_pie = fig_to_image(fig_pie, dpi=130)
    img_pie = Image(buf_pie, width=17*cm, height=14*cm)
    story.append(img_pie)

    # ── PÁGINA 3: Gráfico barras + líneas ────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph('Total reads, Mean Length y N50 por muestra', st_h1))
    story.append(Spacer(1, 4))

    fig_bar = make_barline_chart(df)
    buf_bar = fig_to_image(fig_bar, dpi=140)
    
    # Escalar para ocupar casi toda la página
    img_h = min(26*cm, max(16*cm, n_samples * 0.26 * cm))
    img_bar = Image(buf_bar, width=18*cm, height=img_h)
    story.append(img_bar)

    # Build
    doc.build(story)
    print(f"PDF generado: {output_path}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Uso: python3 {sys.argv[0]} <qc_file.csv> [output.pdf]")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    if len(sys.argv) >= 3:
        out_path = sys.argv[2]
    else:
        base = os.path.splitext(csv_path)[0]
        out_path = base + '_report.pdf'
    
    build_pdf(csv_path, out_path)
