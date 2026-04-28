#!/usr/bin/env python3
"""
build_mag_catalog.py
Construye una tabla consolidada de MAGs integrando:
  - Taxonomia (GTDB-Tk o Sourmash)
  - Calidad (CheckM2)
  - Anotacion (Bakta summary)

Genera: MAGs_Catalog.tsv

Uso:
  python3 scripts/build_mag_catalog.py \
      --taxonomy-dir results/<run>/24_taxonomy_gtdbtk \
      --taxonomy-tool gtdbtk \
      --checkm2-dir results/<run>/23_binqc_checkm2 \
      --bakta-dir results/<run>/31_bakta \
      --output results/<run>/MAGs_Catalog.tsv

  # Con sourmash en vez de GTDB-Tk:
  python3 scripts/build_mag_catalog.py \
      --taxonomy-dir results/<run>/24_taxonomy_sourmash \
      --taxonomy-tool sourmash \
      --checkm2-dir results/<run>/23_binqc_checkm2 \
      --bakta-dir results/<run>/31_bakta \
      --output results/<run>/MAGs_Catalog.tsv
"""

import argparse, os, glob
import pandas as pd

def parse_gtdbtk(tax_dir):
    """Parse GTDB-Tk summary TSV files."""
    rows = []
    for f in sorted(glob.glob(f'{tax_dir}/*/*_taxonomy.tsv')):
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
                        if p.startswith(px):
                            levels[lv] = p[len(px):] or ''
                # Confidence: GTDB-Tk uses ANI or RED values in extra columns
                ani = ''
                for col in df.columns:
                    if 'ani' in col.lower() or 'fastani' in col.lower():
                        ani = str(r.get(col, ''))
                        break
                rows.append({
                    'Sample': sample,
                    'Bin_Name': str(r.iloc[0]).strip(),
                    'Tool': 'GTDB-Tk',
                    'Taxonomic_Confidence': ani or 'pplacer',
                    **levels
                })
        except: pass
    return pd.DataFrame(rows)

def parse_sourmash(tax_dir):
    """Parse Sourmash taxonomy output."""
    rows = []
    for f in sorted(glob.glob(f'{tax_dir}/*/*.taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                bin_name = str(r.get('query_name', r.get('ident', r.iloc[0]))).strip()
                lineage = str(r.get('lineage', ''))
                status = str(r.get('status', ''))
                ani = r.get('query_ani', r.get('ani', ''))
                containment = r.get('f_unique_to_query', r.get('containment', ''))

                levels = {}
                for part in lineage.split(';'):
                    p = part.strip()
                    for px, lv in [('d__','Domain'),('p__','Phylum'),('c__','Class'),
                                   ('o__','Order'),('f__','Family'),('g__','Genus'),('s__','Species')]:
                        if p.startswith(px):
                            levels[lv] = p[len(px):] or ''

                conf = f"ANI={ani}" if ani else f"containment={containment}" if containment else status
                rows.append({
                    'Sample': sample,
                    'Bin_Name': bin_name,
                    'Tool': 'Sourmash',
                    'Taxonomic_Confidence': str(conf),
                    **levels
                })
        except: pass
    return pd.DataFrame(rows)

def parse_checkm2(checkm2_dir):
    rows = []
    for f in sorted(glob.glob(f'{checkm2_dir}/*/quality_report.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                comp = float(r['Completeness']); cont = float(r['Contamination'])
                qc = 'HQ' if comp >= 90 and cont <= 5 else ('MQ' if comp >= 50 and cont <= 10 else 'LQ')
                rows.append({
                    'Sample': sample,
                    'Bin_Name': str(r['Name']).strip(),
                    'Completeness': comp,
                    'Contamination': cont,
                    'Quality': qc
                })
        except: pass
    return pd.DataFrame(rows)

def parse_bakta(bakta_dir):
    rows = []
    for f in sorted(glob.glob(f'{bakta_dir}/*/*/*.txt')):
        bin_name = os.path.basename(os.path.dirname(f))
        sample = os.path.basename(os.path.dirname(os.path.dirname(f)))
        data = {'Sample': sample, 'Bin_Name': bin_name}
        try:
            with open(f) as fh:
                for line in fh:
                    line = line.strip()
                    if ':' not in line: continue
                    key, val = line.split(':', 1)
                    key = key.strip(); val = val.strip()
                    try:
                        if key == 'Length': data['Size_bp'] = int(val)
                        elif key == 'GC': data['GC_pct'] = float(val)
                        elif key == 'tRNAs': data['tRNAs'] = int(val)
                        elif key == 'rRNAs': data['rRNAs'] = int(val)
                        elif key == 'CDSs': data['CDSs'] = int(val)
                        elif key == 'ncRNAs': data['ncRNAs'] = int(val)
                        elif key == 'CRISPR arrays': data['CRISPRs'] = int(val)
                        elif key == 'Count': data['N_contigs'] = int(val)
                        elif key == 'N50': data['N50'] = int(val)
                        elif key.startswith('coding density'): data['Coding_density'] = float(val)
                    except: pass
            rows.append(data)
        except: pass
    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser(description='Build consolidated MAG catalog')
    parser.add_argument('--taxonomy-dir', required=True)
    parser.add_argument('--taxonomy-tool', choices=['gtdbtk', 'sourmash'], default='gtdbtk')
    parser.add_argument('--checkm2-dir', required=True)
    parser.add_argument('--bakta-dir', default='')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # Taxonomy
    if args.taxonomy_tool == 'gtdbtk':
        tax = parse_gtdbtk(args.taxonomy_dir)
    else:
        tax = parse_sourmash(args.taxonomy_dir)
    print(f'[INFO] Taxonomy ({args.taxonomy_tool}): {len(tax)} MAGs')

    # CheckM2
    checkm2 = parse_checkm2(args.checkm2_dir)
    print(f'[INFO] CheckM2: {len(checkm2)} MAGs')

    # Bakta
    bakta = pd.DataFrame()
    if args.bakta_dir and os.path.isdir(args.bakta_dir):
        bakta = parse_bakta(args.bakta_dir)
        print(f'[INFO] Bakta: {len(bakta)} MAGs')

    # Merge
    catalog = checkm2.copy()
    if not tax.empty:
        catalog = catalog.merge(tax, on=['Sample', 'Bin_Name'], how='left')
    if not bakta.empty:
        catalog = catalog.merge(bakta, on=['Sample', 'Bin_Name'], how='left')

    catalog = catalog.fillna('')

    # Order columns
    col_order = ['Sample', 'Bin_Name', 'Completeness', 'Contamination', 'Quality',
                 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
                 'Tool', 'Taxonomic_Confidence',
                 'Size_bp', 'GC_pct', 'N_contigs', 'N50', 'CDSs', 'tRNAs', 'rRNAs',
                 'ncRNAs', 'CRISPRs', 'Coding_density']
    for c in col_order:
        if c not in catalog.columns:
            catalog[c] = ''
    catalog = catalog[col_order]

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    catalog.to_csv(args.output, sep='\t', index=False)
    print(f'[OK] {args.output} ({len(catalog)} MAGs)')

if __name__ == '__main__':
    main()
