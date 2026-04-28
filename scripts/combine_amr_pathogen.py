#!/usr/bin/env python3
"""
combine_amr_pathogen.py
Combina los TSV de integracion AMR-Patogeno-Plasmido de todas las muestras
en un unico informe consolidado.

Uso:
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
    parser.add_argument('--input', nargs='+', required=True, help='TSV files de integracion')
    parser.add_argument('--run-name', required=True, help='Nombre del run')
    parser.add_argument('--output', required=True, help='Output TSV combinado')
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
        print('[WARN] No se encontraron datos de integracion AMR-Patogeno')
        pd.DataFrame().to_csv(args.output, sep='\t', index=False)
        return

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(args.output, sep='\t', index=False)

    # Resumen
    print(f'[INFO] Run: {args.run_name}')
    print(f'[INFO] Muestras: {combined["Sample"].nunique()}')
    print(f'[INFO] MAGs con AMR: {combined["MAG"].nunique()}')
    print(f'[INFO] Total genes AMR: {len(combined)}')
    print(f'[INFO] En plasmido: {len(combined[combined["Location"] == "PLASMID"])}')
    print(f'[INFO] Riesgo CRITICO: {len(combined[combined["Risk_level"] == "CRITICO"])}')
    print(f'[INFO] Riesgo ALTO: {len(combined[combined["Risk_level"] == "ALTO"])}')

    # Top organismos
    print(f'\n[INFO] Top organismos con AMR:')
    for org, count in combined.groupby('Organism')['Gene'].count().sort_values(ascending=False).head(10).items():
        plas = len(combined[(combined['Organism'] == org) & (combined['Location'] == 'PLASMID')])
        print(f'  {org}: {count} genes ({plas} en plasmido)')

    print(f'\n[OK] Reporte combinado: {args.output}')


if __name__ == '__main__':
    main()
