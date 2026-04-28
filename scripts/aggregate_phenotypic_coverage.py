#!/usr/bin/env python3
"""
aggregate_phenotypic_coverage.py
Aggregates minimap2 depth statistics per species from samtools output.
Reports genome coverage breadth at depth thresholds >=1x, >=10x, >=30x.

Input:
  --depth-stats  Per-contig stats from awk (reference, genome_size, bases_1x, bases_10x, bases_30x, depth_sum)
  --idxstats     samtools idxstats output (mapped reads per contig)
  --sample       Sample name

Output: TSV with per-species aggregated coverage metrics.
"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate minimap2 phenotypic coverage per species')
    parser.add_argument('--depth-stats', required=True,
                        help='Per-contig depth stats from awk')
    parser.add_argument('--idxstats', required=True,
                        help='samtools idxstats output')
    parser.add_argument('--sample', required=True,
                        help='Sample name')
    parser.add_argument('--output', required=True,
                        help='Output TSV')
    args = parser.parse_args()

    # Per-contig depth stats (from awk processing of samtools depth -a)
    depth = pd.read_csv(args.depth_stats, sep='\t')

    # Mapped reads per contig (from samtools idxstats)
    idx = pd.read_csv(args.idxstats, sep='\t', header=None,
                       names=['reference', 'length', 'mapped', 'unmapped'])
    idx = idx[idx['reference'] != '*']

    # Extract species name: everything before __ in reference name
    def species_name(ref):
        return ref.split('__')[0].replace('_', ' ')

    depth['species'] = depth['reference'].apply(species_name)
    idx['species'] = idx['reference'].apply(species_name)

    # Aggregate contigs by species (chromosome + plasmids)
    sp = depth.groupby('species').agg({
        'genome_size': 'sum',
        'bases_1x': 'sum',
        'bases_10x': 'sum',
        'bases_30x': 'sum',
        'depth_sum': 'sum'
    }).reset_index()

    reads = idx.groupby('species')['mapped'].sum().reset_index()
    result = sp.merge(reads, on='species', how='left').fillna(0)

    # Compute derived metrics
    result['mean_depth'] = (result['depth_sum'] / result['genome_size']).round(2)
    result['breadth_1x'] = (result['bases_1x'] / result['genome_size'] * 100).round(2)
    result['breadth_10x'] = (result['bases_10x'] / result['genome_size'] * 100).round(2)
    result['breadth_30x'] = (result['bases_30x'] / result['genome_size'] * 100).round(2)
    result['mapped'] = result['mapped'].astype(int)
    result['sample'] = args.sample

    cols = ['sample', 'species', 'genome_size', 'mapped', 'mean_depth',
            'breadth_1x', 'breadth_10x', 'breadth_30x']
    result[cols].to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
