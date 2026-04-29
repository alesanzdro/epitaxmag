#!/usr/bin/env python3
"""
combine_amr_pathogen.py
Combine the AMR-Pathogen-Plasmid integration TSVs from all samples
into a single consolidated report.

Usage:
  python3 combine_amr_pathogen.py \
      --input *.amr_pathogen.tsv \
      --run-name EPIM232 \
      --output combined_amr_pathogen_report.tsv
"""

import argparse
import os
import sys
import glob
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', nargs='+', required=True, help='Integration TSV files')
    parser.add_argument('--run-name', required=True, help='Run name')
    parser.add_argument('--output', required=True, help='Combined output TSV')
    args = parser.parse_args()

    dfs = []
    for f in sorted(args.input):
        if os.path.isfile(f) and os.path.getsize(f) > 0:
            try:
                df = pd.read_csv(f, sep='\t')
                if not df.empty:
                    dfs.append(df)
            except Exception:
                continue

    if not dfs:
        print('[WARN] No AMR-Pathogen integration data found')
        pd.DataFrame().to_csv(args.output, sep='\t', index=False)
        return

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(args.output, sep='\t', index=False)

    # Summary
    print(f'[INFO] Run: {args.run_name}')
    print(f'[INFO] Samples: {combined["Sample"].nunique()}')
    print(f'[INFO] MAGs with AMR: {combined["MAG"].nunique()}')
    print(f'[INFO] Total AMR genes: {len(combined)}')
    print(f'[INFO] On plasmid: {len(combined[combined["Location"] == "PLASMID"])}')
    print(f'[INFO] CRITICAL risk: {len(combined[combined["Risk_level"] == "CRITICAL"])}')
    print(f'[INFO] HIGH risk: {len(combined[combined["Risk_level"] == "HIGH"])}')

    # Top organisms
    print(f'\n[INFO] Top organisms with AMR:')
    for org, count in combined.groupby('Organism')['Gene'].count().sort_values(ascending=False).head(10).items():
        plas = len(combined[(combined['Organism'] == org) & (combined['Location'] == 'PLASMID')])
        print(f'  {org}: {count} genes ({plas} on plasmid)')

    print(f'\n[OK] Combined report: {args.output}')


if __name__ == '__main__':
    main()
