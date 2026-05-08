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
    """Extract genus/species from either a GTDB-Tk string or a Skani/Sourmash
    free-text best hit.

    GTDB-Tk: 'd__Bacteria;p__...;c__...;o__...;f__...;g__Acinetobacter;s__A. baumannii'
    Skani/Sourmash: '<accession> Acinetobacter johnsonii CIP 64.6 contig...'
    """
    levels = {}
    cls = str(classification).strip()
    if not cls:
        return levels

    # GTDB-Tk path: detect by the rank prefixes
    if any(tag in cls for tag in ('d__', 'p__', 'g__', 's__')):
        prefixes = {
            'd__': 'domain', 'p__': 'phylum', 'c__': 'class',
            'o__': 'order', 'f__': 'family', 'g__': 'genus', 's__': 'species'
        }
        for part in cls.split(';'):
            part = part.strip()
            for prefix, level in prefixes.items():
                if part.startswith(prefix):
                    val = part[len(prefix):]
                    levels[level] = val if val else 'Unclassified'
        return levels

    # Skani / Sourmash: '<accession> <Genus species ... description>'
    parts = cls.split(None, 1)
    rest = parts[1] if len(parts) > 1 else ''
    stop_tokens = (' strain ', ' isolate ', ' sp. ', ' subsp. ', ' DSM ',
                   ' ATCC ', ' CIP ', ' NCTC ', ' contig', ' scaffold',
                   ' chromosome', ' complete', ' whole genome', ' MAG:', ',')
    cut = len(rest)
    for tok in stop_tokens:
        i = rest.find(tok)
        if i >= 0 and i < cut:
            cut = i
    binomial = rest[:cut].strip()
    tokens = binomial.split()
    if tokens:
        levels['genus'] = tokens[0]
        if len(tokens) > 1:
            levels['species'] = ' '.join(tokens[1:])
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


# ─── Transparent risk model ────────────────────────────────────
# Stable, ordered rule list. evaluate_risk() walks them top-to-bottom
# and returns the first match plus the rule's human-readable reason.
# The HTML report renders this list verbatim in section "Risk model"
# so users can audit *why* a gene was flagged.

RISK_HIGH_CLASSES = {
    'BETA-LACTAM', 'CARBAPENEM', 'CEPHALOSPORIN',
    'COLISTIN', 'GLYCOPEPTIDE', 'QUINOLONE',
}

RISK_RULES = [
    # (rule_id,    description shown verbatim in the report,
    #  predicate(class_upper, element_type, location, identity, coverage, mag_quality, plasmid_score),
    #  level)
    ('R1',
     'High-risk antibiotic class (carbapenem / β-lactam / colistin / '
     'glycopeptide / quinolone / cephalosporin) on a confirmed plasmid '
     '(geNomad score ≥ 0.7) → highest concern: horizontally transferable '
     'last-resort resistance.',
     lambda cls, et, loc, idn, cov, q, ps: (
         loc == 'PLASMID' and cls in RISK_HIGH_CLASSES and (ps or 0) >= 0.7),
     'CRITICAL'),

    ('R2',
     'High-risk antibiotic class on a chromosome of a high-quality (HQ) MAG '
     'with high sequence support (identity ≥ 95% and coverage ≥ 90%) → '
     'confirmed clinically relevant resistance carrier.',
     lambda cls, et, loc, idn, cov, q, ps: (
         loc == 'CHROMOSOME' and cls in RISK_HIGH_CLASSES
         and q == 'HQ' and (idn or 0) >= 95 and (cov or 0) >= 90),
     'HIGH'),

    ('R3',
     'Any AMR gene on a plasmid (regardless of class) → mobilisable '
     'resistance, monitor.',
     lambda cls, et, loc, idn, cov, q, ps: loc == 'PLASMID',
     'HIGH'),

    ('R4',
     'High-risk class on a chromosome but evidence is weaker '
     '(MQ/LQ MAG or low identity/coverage) → flagged but needs '
     'confirmation.',
     lambda cls, et, loc, idn, cov, q, ps: (
         cls in RISK_HIGH_CLASSES),
     'MEDIUM'),

    ('R5',
     'Default: any other AMR gene call.',
     lambda cls, et, loc, idn, cov, q, ps: True,
     'LOW'),
]


def evaluate_risk(amr_class, element_type, location,
                  identity=None, coverage=None,
                  mag_quality=None, plasmid_score=None):
    """Apply the ordered RISK_RULES and return (level, rule_id, reason)."""
    cls = (amr_class or '').upper()
    for rule_id, reason, pred, level in RISK_RULES:
        try:
            if pred(cls, element_type, location, identity, coverage,
                    mag_quality, plasmid_score):
                return level, rule_id, reason
        except Exception:
            continue
    return 'LOW', 'R5', RISK_RULES[-1][1]


def classify_risk(amr_class, element_type, location):
    """Backwards-compatible wrapper used by older callers."""
    level, _, _ = evaluate_risk(amr_class, element_type, location)
    return level


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

        # Transparent risk evaluation (level + rule_id + reason)
        risk, risk_rule, risk_reason = evaluate_risk(
            hit['amr_class'], hit['element_type'], location,
            identity=hit['pct_identity'], coverage=hit['pct_coverage'],
            mag_quality=qc_label, plasmid_score=plas_score)

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
            'Risk_rule': risk_rule,
            'Risk_reason': risk_reason,
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
