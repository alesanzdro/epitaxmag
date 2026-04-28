#!/usr/bin/env python3
"""
build_phenotypic_db.py
Downloads type-strain genomes from NCBI and creates a combined FASTA
for phenotypic verification with minimap2.

Input:  TSV with columns 'organism' and 'accession'
Output: Combined FASTA with headers >Species_name__Accession

Uso:
    python3 scripts/build_phenotypic_db.py \
        --targets assets/phenotypic_targets/target_organisms.tsv \
        --output phenotypic_targets.fa
"""

import argparse, csv, os, sys, tempfile, shutil, zipfile, urllib.request


def download_genome(accession, species_tag, tmpdir):
    """Download genome FASTA from NCBI datasets API and rewrite headers."""
    zip_path = os.path.join(tmpdir, f'{species_tag}.zip')

    url = (f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'
           f'{accession}/download?include_annotation_type=GENOME_FASTA')

    try:
        urllib.request.urlretrieve(url, zip_path)
    except Exception as e:
        print(f'  Download error: {e}', file=sys.stderr)
        return None

    if not os.path.isfile(zip_path) or os.path.getsize(zip_path) == 0:
        return None

    # Extract .fna files from zip (NCBI datasets format)
    extract_dir = os.path.join(tmpdir, f'extract_{species_tag}')
    os.makedirs(extract_dir, exist_ok=True)

    try:
        with zipfile.ZipFile(zip_path) as zf:
            fna_members = [m for m in zf.namelist() if m.endswith('.fna')]
            for m in fna_members:
                zf.extract(m, extract_dir)
    except Exception as e:
        print(f'  Extraction error: {e}', file=sys.stderr)
        return None

    # Find all extracted .fna files
    fna_files = []
    for root, _, files in os.walk(extract_dir):
        for f in sorted(files):
            if f.endswith('.fna'):
                fna_files.append(os.path.join(root, f))

    if not fna_files:
        return None

    # Read sequences and rewrite headers as Species_name__ContigAccession
    lines = []
    for fna_path in fna_files:
        with open(fna_path) as f:
            for line in f:
                if line.startswith('>'):
                    contig_id = line[1:].split()[0]
                    lines.append(f'>{species_tag}__{contig_id}\n')
                else:
                    lines.append(line)

    return ''.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Build phenotypic verification FASTA from NCBI type-strain genomes')
    parser.add_argument('--targets', required=True,
                        help='TSV with columns: organism, accession')
    parser.add_argument('--output', required=True,
                        help='Output combined FASTA')
    args = parser.parse_args()

    # Read targets TSV
    targets = []
    with open(args.targets) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            targets.append((row['organism'], row['accession']))

    print(f'[INFO] Downloading {len(targets)} type-strain genomes...', file=sys.stderr)

    ok = 0
    with open(args.output, 'w') as out:
        for organism, accession in targets:
            species_tag = organism.replace(' ', '_')
            print(f'  {organism} ({accession})...', end=' ', file=sys.stderr, flush=True)

            with tempfile.TemporaryDirectory() as tmpdir:
                content = download_genome(accession, species_tag, tmpdir)

            if content:
                out.write(content)
                print('OK', file=sys.stderr)
                ok += 1
            else:
                print('FAILED', file=sys.stderr)

    n_contigs = sum(1 for line in open(args.output) if line.startswith('>'))
    print(f'[INFO] Database: {n_contigs} contigs from {ok}/{len(targets)} species',
          file=sys.stderr)

    if ok == 0:
        print('[ERROR] No genomes downloaded. Check network/proxy.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
