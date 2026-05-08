#!/usr/bin/env python3
"""
Comprehensive analysis of Dorado FASTQ files.
Reports sequence metrics, ASCII quality, Dorado quality (qs:f:) and a
custom global quality score (q_custom_mean).
"""

import os
import gzip
import argparse
from concurrent.futures import ProcessPoolExecutor
import time
import pandas as pd
import numpy as np
import statistics
import math

def phred_to_prob_error(q):
    return 10 ** (-q / 10)

def prob_error_to_phred(p):
    return -10 * math.log10(p) if p > 0 else float('inf')

def calculate_q_custom_mean(ascii_quals):
    """
    Compute the global per-read quality (q_custom_mean):
      1. Convert each Q to its error probability.
      2. Take the arithmetic mean of those probabilities.
      3. Convert the mean back to Phred.
    """
    probs = [phred_to_prob_error(q) for q in ascii_quals]
    p_avg = sum(probs) / len(probs)
    return prob_error_to_phred(p_avg)

def analyze_fastq_complete(file_path, tagdelete):
    print(f"Procesando archivo: {os.path.basename(file_path)}")

    # contenedores
    simplex_dorado_quals = []
    duplex_dorado_quals = []
    all_lengths = []
    all_ascii_quals = []
    all_custom_quals = []
    total_gc_count = 0
    total_bases = 0
    read_count = 0

    with gzip.open(file_path, 'rt') as f:
        lines = []
        for line in f:
            lines.append(line.strip())
            if len(lines) == 4:
                header, sequence, _, quality_str = lines
                read_count += 1

                # q-header Dorado
                dorado_q = None
                for part in header.split():  # Cambio: usar split() sin argumentos para separar por espacios
                    if part.startswith('qs:f:'):
                        try:
                            dorado_q = float(part.split(':')[2])
                        except:
                            pass

                # longitud y GC
                L = len(sequence)
                all_lengths.append(L)
                total_bases += L
                total_gc_count += sequence.count('G') + sequence.count('C')

                # ASCII quals
                ascii_quals = [ord(c)-33 for c in quality_str]
                avg_ascii = sum(ascii_quals) / len(ascii_quals) if ascii_quals else 0
                all_ascii_quals.append(avg_ascii)

                # q_custom_mean
                q_custom = calculate_q_custom_mean(ascii_quals) if ascii_quals else 0
                all_custom_quals.append(q_custom)

                # almacenar q-header
                if dorado_q is not None:
                    if ';' in header:
                        duplex_dorado_quals.append(dorado_q)
                    else:
                        simplex_dorado_quals.append(dorado_q)

                lines = []

    # basic metrics
    all_lengths.sort()
    total_reads = read_count
    min_len = all_lengths[0] if all_lengths else 0
    max_len = all_lengths[-1] if all_lengths else 0
    mean_len = statistics.mean(all_lengths) if all_lengths else 0
    med_len = statistics.median(all_lengths) if all_lengths else 0

    # N50
    n50 = 0
    if all_lengths:
        half = sum(all_lengths)/2
        cum = 0
        for L in sorted(all_lengths, reverse=True):
            cum += L
            if cum >= half:
                n50 = L
                break

    # mean values
    all_dorado_quals = simplex_dorado_quals + duplex_dorado_quals
    mean_dorado = statistics.mean(all_dorado_quals) if all_dorado_quals else 0
    mean_ascii = statistics.mean(all_ascii_quals) if all_ascii_quals else 0
    mean_custom = statistics.mean(all_custom_quals) if all_custom_quals else 0

    # build output
    sample = os.path.basename(file_path).replace(tagdelete, '').replace('.fastq.gz','')
    return sample, {
        'total_reads': total_reads,
        'min_length': min_len,
        'max_length': max_len,
        'mean_length': mean_len,
        'median_length': med_len,
        'n50': n50,
        'gc_percentage': (total_gc_count/total_bases*100) if total_bases else 0,
        'dorado_mean_q': mean_dorado,
        'ascii_mean_q': mean_ascii,
        'custom_mean_q': mean_custom
    }

def main():
    parser = argparse.ArgumentParser(description="Full FASTQ analysis with q_custom_mean")
    parser.add_argument("dir", help="Directory containing .fastq.gz files")
    parser.add_argument("-t","--threads", type=int, default=4, help="Worker threads")
    parser.add_argument("--tagdelete", default="", help="Substring to strip from sample names")
    args = parser.parse_args()

    files = [os.path.join(args.dir,f) for f in os.listdir(args.dir) if f.endswith('.fastq.gz')]
    results = []
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        for res in ex.map(analyze_fastq_complete, files, [args.tagdelete]*len(files)):
            if res: results.append(res)

    # DataFrame
    df = pd.DataFrame([ {'Sample':s, **m} for s,m in results ])
    df = df[['Sample','total_reads','min_length','max_length','mean_length','median_length','n50',
             'gc_percentage','dorado_mean_q','ascii_mean_q','custom_mean_q']]
    out_csv = os.path.join(args.dir, 'analysis_with_q_custom.csv')
    df.to_csv(out_csv, index=False, sep=';', decimal=',')

    print(f"Guardado: {out_csv}")

if __name__ == "__main__":
    main()
