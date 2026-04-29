#!/usr/bin/env python3
"""
phenotypic_check.py
Cross conventional-microbiology phenotypic results with metagenomic
evidence from multiple sources to verify target organisms.

Evidence cascade:
  Tier 1: MAG with GTDB-Tk match (completeness >=50%)
  Tier 2: Sylph hit (ANI >=95, abundance >=0.01%)
  Tier 3: Kraken2/Bracken (reads >=50)
  Tier 4: KMA type-strain (depth >=1x, identity >=90%)
  Tier 0: not detected

Usage:
  python3 scripts/phenotypic_check.py \
      --results-dir results/260226_EPIM232 \
      --phenotypic assets/phenotypic_targets/260226_EPIM232.tsv \
      --kma-phenotypic results/260226_EPIM232/33_phenotypic_kma \
      --output results/260226_EPIM232/32_reports/phenotypic_summary.tsv
"""

import argparse, os, glob
import pandas as pd
import numpy as np

# Thresholds
SYLPH_MIN_ANI = 95.0
SYLPH_MIN_ABUND = 0.01  # %
BRACKEN_MIN_READS = 50
KMA_MIN_DEPTH = 1.0
KMA_MIN_IDENTITY = 90.0
CHECKM2_MIN_COMP = 50.0


def load_phenotypic(path):
    if not os.path.isfile(path): return pd.DataFrame()
    return pd.read_csv(path, sep='\t')


def load_sylph(results_dir):
    path = f'{results_dir}/10_tax_sylph/sylph_profile_all.tsv'
    if not os.path.isfile(path): return pd.DataFrame()
    df = pd.read_csv(path, sep='\t')
    df['Sample'] = df['Sample_file'].apply(
        lambda s: os.path.basename(str(s)).replace('.filtered.fastq.gz','').replace('.clean.fastq.gz',''))
    df['Tax_abund'] = pd.to_numeric(df['Taxonomic_abundance'], errors='coerce').fillna(0)
    df['Adj_ANI'] = pd.to_numeric(df['Adjusted_ANI'], errors='coerce').fillna(0)
    df['Eff_cov'] = pd.to_numeric(df['Eff_cov'], errors='coerce').fillna(0)
    df['Contig_name'] = df['Contig_name'].astype(str)
    return df


def load_bracken(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/08_tax_bracken/*.bracken.S.txt')):
        sample = os.path.basename(f).replace('.bracken.S.txt', '')
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample,
                    'Name': str(r.get('name', '')).strip(),
                    'Reads': int(r.get('new_est_reads', r.get('kraken_assigned_reads', 0))),
                    'Fraction': float(r.get('fraction_total_reads', 0))
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_gtdbtk(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/24_taxonomy_gtdbtk/*/*_taxonomy.tsv') +
                     glob.glob(f'{results_dir}/24_taxonomy_sourmash/*/*_taxonomy.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                cls = str(r.iloc[1]) if len(r) > 1 else ''
                rows.append({'Sample': sample, 'Bin': str(r.iloc[0]).strip(), 'Taxonomy': cls})
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_checkm2(results_dir):
    rows = []
    for f in sorted(glob.glob(f'{results_dir}/23_binqc_checkm2/*/quality_report.tsv')):
        sample = os.path.basename(os.path.dirname(f))
        try:
            df = pd.read_csv(f, sep='\t')
            for _, r in df.iterrows():
                rows.append({
                    'Sample': sample, 'Bin': str(r['Name']).strip(),
                    'Completeness': float(r['Completeness']),
                    'Contamination': float(r['Contamination'])
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_kma_phenotypic(kma_dir):
    rows = []
    if not kma_dir or not os.path.isdir(kma_dir): return pd.DataFrame()
    for f in sorted(glob.glob(f'{kma_dir}/*.res')):
        sample = os.path.basename(f).replace('_phenotypic.res', '').replace('.res', '')
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty: continue
            for _, r in df.iterrows():
                template = str(r.iloc[0]).strip()
                if not template or template.startswith('#'): continue
                # Template format: "Species_name contig_info"
                species = template.split()[0].replace('_', ' ') if template else ''
                rows.append({
                    'Sample': sample,
                    'Template': template,
                    'Species_kma': species,
                    'Depth': float(r.get('Depth', 0)),
                    'Identity': float(r.get('Template_Identity', 0)),
                    'Coverage': float(r.get('Template_Coverage', 0))
                })
        except: pass
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def search_sylph(sylph_df, sample, organism):
    """Search Sylph for a target organism (genus or species match)."""
    if sylph_df.empty: return None
    genus = organism.split()[0]
    species = organism.lower()

    mask = (sylph_df['Sample'] == sample) & (sylph_df['Contig_name'].str.contains(genus, case=False, na=False))
    hits = sylph_df[mask].copy()

    if hits.empty: return None

    # Try species-level match first
    sp_hits = hits[hits['Contig_name'].str.lower().str.contains(species, na=False)]
    if not sp_hits.empty:
        best = sp_hits.loc[sp_hits['Tax_abund'].idxmax()]
    else:
        best = hits.loc[hits['Tax_abund'].idxmax()]

    return {
        'Sylph_abundance': round(best['Tax_abund'], 4),
        'Sylph_ANI': round(best['Adj_ANI'], 2),
        'Sylph_coverage': round(best['Eff_cov'], 3),
        'Sylph_match': 'species' if not sp_hits.empty else 'genus'
    }


def search_bracken(bracken_df, sample, organism):
    """Search Bracken for a target organism."""
    if bracken_df.empty: return None
    genus = organism.split()[0]
    species = organism.lower()

    mask = (bracken_df['Sample'] == sample) & (bracken_df['Name'].str.contains(genus, case=False, na=False))
    hits = bracken_df[mask]

    if hits.empty: return None

    sp_hits = hits[hits['Name'].str.lower().str.contains(species, na=False)]
    if not sp_hits.empty:
        best = sp_hits.loc[sp_hits['Reads'].idxmax()]
    else:
        best = hits.loc[hits['Reads'].idxmax()]

    return {
        'Bracken_reads': int(best['Reads']),
        'Bracken_abundance': round(best['Fraction'] * 100, 4),
        'Bracken_name': best['Name'],
        'Bracken_match': 'species' if not sp_hits.empty else 'genus'
    }


def search_mag(gtdbtk_df, checkm2_df, sample, organism):
    """Search for a MAG matching the target organism."""
    if gtdbtk_df.empty: return None
    genus = organism.split()[0]

    mask = (gtdbtk_df['Sample'] == sample) & (gtdbtk_df['Taxonomy'].str.contains(genus, case=False, na=False))
    hits = gtdbtk_df[mask]

    if hits.empty: return None

    # Check completeness
    results = []
    for _, h in hits.iterrows():
        bin_name = h['Bin']
        qc = checkm2_df[(checkm2_df['Sample'] == sample) & (checkm2_df['Bin'] == bin_name)]
        comp = float(qc['Completeness'].iloc[0]) if not qc.empty else 0
        cont = float(qc['Contamination'].iloc[0]) if not qc.empty else 0
        if comp >= CHECKM2_MIN_COMP:
            results.append({
                'MAG_name': bin_name,
                'MAG_completeness': round(comp, 1),
                'MAG_contamination': round(cont, 1),
                'MAG_taxonomy': h['Taxonomy']
            })

    if results:
        return max(results, key=lambda x: x['MAG_completeness'])
    return None


def search_kma(kma_df, sample, organism):
    """Search KMA type-strain results."""
    if kma_df.empty: return None
    genus = organism.split()[0]
    species = organism.replace(' ', '_').lower()

    mask = (kma_df['Sample'] == sample) & (
        kma_df['Template'].str.lower().str.contains(genus.lower(), na=False))
    hits = kma_df[mask]

    if hits.empty: return None

    # Best hit by depth
    best = hits.loc[hits['Depth'].idxmax()]

    return {
        'KMA_depth': round(best['Depth'], 2),
        'KMA_identity': round(best['Identity'], 1),
        'KMA_coverage': round(best['Coverage'], 1),
        'KMA_template': best['Template'][:50]
    }


def determine_tier(mag_result, sylph_result, bracken_result, kma_result):
    """Determine best evidence tier."""
    if mag_result and mag_result['MAG_completeness'] >= CHECKM2_MIN_COMP:
        return 1, 'MAG Confirmado'
    if sylph_result and sylph_result['Sylph_ANI'] >= SYLPH_MIN_ANI and sylph_result['Sylph_abundance'] >= SYLPH_MIN_ABUND:
        return 2, 'Sylph Detectado'
    if bracken_result and bracken_result['Bracken_reads'] >= BRACKEN_MIN_READS:
        return 3, 'Bracken Trazas'
    if kma_result and kma_result['KMA_depth'] >= KMA_MIN_DEPTH and kma_result['KMA_identity'] >= KMA_MIN_IDENTITY:
        return 4, 'KMA Read-level'
    # Check with relaxed thresholds
    if sylph_result and sylph_result['Sylph_abundance'] > 0:
        return 3, 'Sylph Trazas'
    if bracken_result and bracken_result['Bracken_reads'] > 0:
        return 3, 'Bracken Trazas (bajo)'
    if kma_result and kma_result['KMA_depth'] > 0:
        return 4, 'KMA Trazas'
    return 0, 'No detectado'


def determine_concordance(pheno_detected, tier):
    if pheno_detected == 'unknown':
        return 'unknown'
    pheno = pheno_detected == 'yes'
    seq = tier > 0
    if pheno and seq: return 'agree'
    if pheno and not seq: return 'only_phenotypic'
    if not pheno and seq: return 'only_sequencing'
    return 'both_negative'


def main():
    parser = argparse.ArgumentParser(description='Phenotypic vs sequencing verification')
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--phenotypic', required=True)
    parser.add_argument('--kma-phenotypic', default='')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    rd = args.results_dir
    print('[INFO] Loading data sources...')

    pheno = load_phenotypic(args.phenotypic)
    print(f'  Phenotypic: {len(pheno)} entries')

    sylph = load_sylph(rd)
    print(f'  Sylph: {len(sylph)} hits')

    bracken = load_bracken(rd)
    print(f'  Bracken: {len(bracken)} entries')

    gtdbtk = load_gtdbtk(rd)
    print(f'  GTDB-Tk/Sourmash: {len(gtdbtk)} MAGs')

    checkm2 = load_checkm2(rd)
    print(f'  CheckM2: {len(checkm2)} MAGs')

    kma = load_kma_phenotypic(args.kma_phenotypic)
    print(f'  KMA type-strain: {len(kma)} hits')

    # Process each phenotypic entry
    rows = []
    for _, p in pheno.iterrows():
        sample = p['sample_id']
        organism = p['organism']

        mag_result = search_mag(gtdbtk, checkm2, sample, organism)
        sylph_result = search_sylph(sylph, sample, organism)
        bracken_result = search_bracken(bracken, sample, organism)
        kma_result = search_kma(kma, sample, organism)

        tier, evidence = determine_tier(mag_result, sylph_result, bracken_result, kma_result)
        concordance = determine_concordance(p['detected'], tier)

        row = {
            'Sample': sample,
            'Organism': organism,
            'Phenotypic': p['detected'],
            'Method': p['method'],
            'Notes': p.get('notes', ''),
            'Seq_tier': tier,
            'Seq_evidence': evidence,
            'Concordance': concordance,
        }

        # Add details from each source
        if mag_result:
            row.update({k: v for k, v in mag_result.items()})
        else:
            row.update({'MAG_name': '', 'MAG_completeness': '', 'MAG_contamination': '', 'MAG_taxonomy': ''})

        if sylph_result:
            row.update(sylph_result)
        else:
            row.update({'Sylph_abundance': '', 'Sylph_ANI': '', 'Sylph_coverage': '', 'Sylph_match': ''})

        if bracken_result:
            row.update(bracken_result)
        else:
            row.update({'Bracken_reads': '', 'Bracken_abundance': '', 'Bracken_name': '', 'Bracken_match': ''})

        if kma_result:
            row.update(kma_result)
        else:
            row.update({'KMA_depth': '', 'KMA_identity': '', 'KMA_coverage': '', 'KMA_template': ''})

        rows.append(row)

    result = pd.DataFrame(rows)

    # Column order
    cols = ['Sample', 'Organism', 'Phenotypic', 'Method', 'Seq_tier', 'Seq_evidence', 'Concordance',
            'MAG_name', 'MAG_completeness', 'MAG_contamination',
            'Sylph_abundance', 'Sylph_ANI', 'Sylph_coverage',
            'Bracken_reads', 'Bracken_abundance',
            'KMA_depth', 'KMA_identity', 'KMA_coverage',
            'Notes']
    for c in cols:
        if c not in result.columns: result[c] = ''
    result = result[cols]

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    result.to_csv(args.output, sep='\t', index=False)

    # Summary
    print(f'\n[RESULTS] {len(result)} organism x sample combinations')
    for conc in ['agree', 'only_phenotypic', 'only_sequencing', 'both_negative', 'unknown']:
        n = len(result[result['Concordance'] == conc])
        if n > 0: print(f'  {conc}: {n}')

    print(f'\n[OK] {args.output}')


if __name__ == '__main__':
    main()
