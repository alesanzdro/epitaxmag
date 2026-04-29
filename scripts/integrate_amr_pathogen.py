#!/usr/bin/env python3
"""
integrate_amr_pathogen.py
AMR-Pathogen-Plasmid integration for EpiTaxMAG.

Cross-references results from:
  - GTDB-Tk (MAG taxonomy)
  - AMRFinderPlus (AMR genes per bin)
  - geNomad (plasmid/viral contigs)
  - CheckM2 (bin quality)

Outputs:
  - Integrated TSV with: organism + AMR gene + location (chromosome/plasmid)
  - Risk summary per organism

Usage:
  python3 integrate_amr_pathogen.py \
      --sample SAMPLE_NAME \
      --taxonomy sample_taxonomy.tsv \
      --amr sample_amr_combined.tsv \
      --checkm2 quality_report.tsv \
      --plasmids sample_plasmid_summary.tsv \
      --output sample_amr_pathogen.tsv
"""

import argparse
import os
import sys
import glob
import pandas as pd
from datetime import datetime


def parse_taxonomy(path):
    """Parse GTDB-Tk summary TSV -> dict[bin_name] = taxonomy_string."""
    tax = {}
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return tax
    try:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            return tax
        # GTDB-Tk output: user_genome, classification, ...
        col_genome = df.columns[0]
        col_class = 'classification' if 'classification' in df.columns else df.columns[1]
        for _, row in df.iterrows():
            bin_name = str(row[col_genome]).strip()
            classification = str(row[col_class]).strip()
            tax[bin_name] = classification
    except Exception:
        pass
    return tax


def parse_taxonomy_levels(classification):
    """Extract taxonomic levels from a GTDB string.
    Input: d__Bacteria;p__Proteobacteria;c__...;g__Acinetobacter;s__A. baumannii
    Output: dict with domain, phylum, class, order, family, genus, species
    """
    levels = {}
    prefixes = {
        'd__': 'domain', 'p__': 'phylum', 'c__': 'class',
        'o__': 'order', 'f__': 'family', 'g__': 'genus', 's__': 'species'
    }
    for part in str(classification).split(';'):
        part = part.strip()
        for prefix, level in prefixes.items():
            if part.startswith(prefix):
                val = part[len(prefix):]
                levels[level] = val if val else 'Unclassified'
    return levels


def parse_checkm2(path):
    """Parse CheckM2 quality_report.tsv -> dict[bin_name] = (completeness, contamination)."""
    qc = {}
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return qc
    try:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            return qc
        for _, row in df.iterrows():
            name = str(row['Name']).strip()
            comp = float(row['Completeness'])
            cont = float(row['Contamination'])
            qc[name] = (comp, cont)
    except Exception:
        pass
    return qc


def parse_amr(path):
    """Parse AMRFinderPlus results -> list of dicts.
    Accepts a directory (per_bin/*.tsv) or a combined file.
    """
    hits = []

    # If it's a directory, read individual files
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, '*_amr.tsv')))
    elif os.path.isfile(path):
        # If the combined file has a corrupt header, try per_bin/
        parent = os.path.dirname(path)
        per_bin = os.path.join(parent, 'per_bin')
        if os.path.isdir(per_bin):
            files = sorted(glob.glob(os.path.join(per_bin, '*_amr.tsv')))
        else:
            files = [path]
    else:
        return hits

    for fpath in files:
        try:
            df = pd.read_csv(fpath, sep='\t', comment='=')
            if df.empty:
                continue
            for _, row in df.iterrows():
                hit = {
                    'bin_name': str(row.get('Name', '')).strip(),
                    'contig_id': str(row.get('Contig id', '')).strip(),
                    'gene_symbol': str(row.get('Element symbol', row.get('Gene symbol', ''))).strip(),
                    'sequence_name': str(row.get('Element name', row.get('Sequence name', ''))).strip(),
                    'scope': str(row.get('Scope', '')).strip(),
                    'element_type': str(row.get('Type', row.get('Element type', ''))).strip(),
                    'element_subtype': str(row.get('Subtype', row.get('Element subtype', ''))).strip(),
                    'amr_class': str(row.get('Class', '')).strip(),
                    'amr_subclass': str(row.get('Subclass', '')).strip(),
                    'method': str(row.get('Method', '')).strip(),
                    'pct_identity': float(row.get('% Identity to reference', row.get('% Identity to reference sequence', 0))),
                    'pct_coverage': float(row.get('% Coverage of reference', row.get('% Coverage of reference sequence', 0))),
                }
                hits.append(hit)
        except Exception:
            continue
    return hits


def parse_plasmids(path):
    """Parse geNomad plasmid_summary.tsv -> dict[contig_name] = plasmid_info."""
    plasmids = {}
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return plasmids
    try:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            return plasmids
        for _, row in df.iterrows():
            seq_name = str(row.iloc[0]).strip()
            # geNomad appends |provirus or similar, clean it up
            base_name = seq_name.split('|')[0]
            score = float(row.get('plasmid_score', row.get('score', 0)))
            n_genes = int(row.get('n_genes', 0)) if 'n_genes' in row else 0
            plasmids[base_name] = {
                'plasmid_score': score,
                'n_genes_plasmid': n_genes,
            }
    except Exception:
        pass
    return plasmids


def classify_risk(amr_class, element_type, location):
    """Classify risk level of an AMR gene."""
    high_risk_classes = {
        'BETA-LACTAM', 'CARBAPENEM', 'CEPHALOSPORIN',
        'COLISTIN', 'GLYCOPEPTIDE', 'QUINOLONE'
    }
    amr_upper = amr_class.upper() if amr_class else ''

    if location == 'PLASMID' and amr_upper in high_risk_classes:
        return 'CRITICAL'
    elif location == 'PLASMID':
        return 'HIGH'
    elif amr_upper in high_risk_classes:
        return 'MEDIUM'
    else:
        return 'LOW'


def quality_label(comp, cont):
    """MIMAG quality label."""
    if comp >= 90 and cont <= 5:
        return 'HQ'
    elif comp >= 50 and cont <= 10:
        return 'MQ'
    else:
        return 'LQ'


def integrate(sample, taxonomy, checkm2, amr_hits, plasmids):
    """Integrate all sources into a unified table."""
    rows = []

    for hit in amr_hits:
        bin_name = hit['bin_name']
        contig = hit['contig_id']

        # Taxonomy
        tax_string = taxonomy.get(bin_name, '')
        levels = parse_taxonomy_levels(tax_string)
        genus = levels.get('genus', 'Unknown')
        species = levels.get('species', '')
        organism = f"{genus} {species}".strip() if species else genus

        # Bin quality
        comp, cont = checkm2.get(bin_name, (0, 0))
        qc_label = quality_label(comp, cont)

        # Location: plasmid or chromosome?
        if contig in plasmids:
            location = 'PLASMID'
            plas_score = plasmids[contig]['plasmid_score']
        else:
            location = 'CHROMOSOME'
            plas_score = 0.0

        # Risk
        risk = classify_risk(hit['amr_class'], hit['element_type'], location)

        rows.append({
            'Sample': sample,
            'MAG': bin_name,
            'Organism': organism,
            'Taxonomy': tax_string,
            'MAG_quality': qc_label,
            'Completeness': round(comp, 1),
            'Contamination': round(cont, 1),
            'Gene': hit['gene_symbol'],
            'Gene_description': hit['sequence_name'],
            'AMR_class': hit['amr_class'],
            'AMR_subclass': hit['amr_subclass'],
            'Element_type': hit['element_type'],
            'Method': hit['method'],
            'Identity_pct': round(hit['pct_identity'], 1),
            'Coverage_pct': round(hit['pct_coverage'], 1),
            'Contig': contig,
            'Location': location,
            'Plasmid_score': round(plas_score, 3),
            'Risk_level': risk,
        })

    return pd.DataFrame(rows)


def generate_summary(df):
    """Generate a text summary from the integrated table."""
    if df.empty:
        return "No resistance genes detected in the MAGs of this sample."

    lines = []
    sample = df['Sample'].iloc[0]
    n_genes = len(df)
    n_mags = df['MAG'].nunique()
    n_plasmid = len(df[df['Location'] == 'PLASMID'])
    n_critical = len(df[df['Risk_level'] == 'CRITICAL'])
    n_high = len(df[df['Risk_level'] == 'HIGH'])

    lines.append(f"Sample: {sample}")
    lines.append(f"Total AMR genes: {n_genes} in {n_mags} MAGs")
    lines.append(f"On plasmids: {n_plasmid} ({n_plasmid/n_genes*100:.0f}%)")

    if n_critical > 0:
        lines.append(f"ALERT: {n_critical} CRITICAL AMR genes (key resistance on plasmid)")
        critical = df[df['Risk_level'] == 'CRITICAL']
        for _, row in critical.iterrows():
            lines.append(f"  - {row['Gene']} ({row['AMR_class']}) in {row['Organism']} [plasmid, score={row['Plasmid_score']:.2f}]")

    if n_high > 0:
        lines.append(f"HIGH risk: {n_high} AMR genes on plasmids (transferable)")

    # Top organisms
    lines.append("")
    lines.append("Organisms with AMR:")
    org_counts = df.groupby('Organism')['Gene'].count().sort_values(ascending=False)
    for org, count in org_counts.items():
        sub = df[df['Organism'] == org]
        classes = ', '.join(sorted(sub['AMR_class'].unique()))
        plas = len(sub[sub['Location'] == 'PLASMID'])
        lines.append(f"  {org}: {count} genes ({classes}) [{plas} on plasmid]")

    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--taxonomy', required=True, help='GTDB-Tk taxonomy TSV')
    parser.add_argument('--amr', required=True, help='AMRFinderPlus combined TSV')
    parser.add_argument('--checkm2', required=True, help='CheckM2 quality_report.tsv')
    parser.add_argument('--plasmids', default='', help='geNomad plasmid_summary.tsv')
    parser.add_argument('--output', required=True, help='Output integrated TSV')
    parser.add_argument('--summary', default='', help='Output summary TXT')
    args = parser.parse_args()

    print(f'[INFO] Integrating AMR-Pathogen-Plasmid for {args.sample}')

    taxonomy = parse_taxonomy(args.taxonomy)
    print(f'[INFO] GTDB-Tk: {len(taxonomy)} MAGs classified')

    checkm2 = parse_checkm2(args.checkm2)
    print(f'[INFO] CheckM2: {len(checkm2)} bins evaluated')

    amr_hits = parse_amr(args.amr)
    print(f'[INFO] AMRFinderPlus: {len(amr_hits)} AMR genes')

    plasmids = parse_plasmids(args.plasmids)
    print(f'[INFO] geNomad: {len(plasmids)} plasmid contigs')

    df = integrate(args.sample, taxonomy, checkm2, amr_hits, plasmids)

    # Save TSV
    df.to_csv(args.output, sep='\t', index=False)
    print(f'[OK] Integrated table: {args.output} ({len(df)} rows)')

    # Save summary
    summary = generate_summary(df)
    summary_path = args.summary or args.output.replace('.tsv', '_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(summary)
    print(f'[OK] Summary: {summary_path}')
    print(summary)


if __name__ == '__main__':
    main()
