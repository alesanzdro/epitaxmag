#!/usr/bin/env nextflow

/*
 * SETUP_DB — Automatic download of pipeline databases.
 *
 * Each process uses `storeDir` so files persist across runs and across
 * pipeline restarts. Skipped automatically if the marker file already
 * exists. Use `--skip_db_setup` to bypass entirely.
 *
 * Layout under params.db_root:
 *   k2_pluspf_20251015/        — Kraken2 PlusPF
 *   kaiju/                     — Kaiju RefSeq (FMI + nodes/names dump)
 *   sylph/                     — Sylph GTDB + Fungi + Viral sketches
 *   kma_resfinder/resfinder_db/ — ResFinder concatenated FASTA
 *   fastqscreen/               — FastQ Screen panel + config
 *   checkm2/CheckM2_database/  — CheckM2 DIAMOND DB
 *   gtdbtk_r232/release232/    — GTDB-Tk r232 reference data
 *   genomad/genomad_db/        — geNomad v1.9 database
 *   amrfinderplus/latest/      — AMRFinderPlus DB (downloaded by tool)
 *   bakta/db/                  — Bakta DB (downloaded by tool)
 *   sourmash/                  — Sourmash GTDB sketches (only if --skip_gtdbtk)
 */


// ──────────────────────────────────────────────────────────────────────
// TAX-phase databases
// ──────────────────────────────────────────────────────────────────────

process DOWNLOAD_KRAKEN2 {
    storeDir "${params.db_root}/k2_pluspf_20251015"

    output:
    path("hash.k2d"), emit: done

    when:
    !file("${params.db_root}/k2_pluspf_20251015/hash.k2d").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading Kraken2 PlusPF (~78 GB compressed)..."
    wget -c --tries=5 --timeout=120 \\
        https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz
    tar xzf k2_pluspf_20251015.tar.gz
    rm -f k2_pluspf_20251015.tar.gz
    """
}

process DOWNLOAD_KAIJU {
    storeDir "${params.db_root}/kaiju"

    output:
    path("kaiju_db_refseq_ref.fmi"), emit: done

    when:
    !file("${params.db_root}/kaiju/kaiju_db_refseq_ref.fmi").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading Kaiju RefSeq ref index 2024-08-14 (~36 GB compressed)..."
    wget -c --tries=5 --timeout=120 \\
        https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz
    tar xzf kaiju_db_refseq_ref_2024-08-14.tgz
    rm -f kaiju_db_refseq_ref_2024-08-14.tgz
    # Sanity check: the expected files must exist
    test -f kaiju_db_refseq_ref.fmi
    test -f nodes.dmp
    test -f names.dmp
    """
}

process DOWNLOAD_SYLPH {
    storeDir "${params.db_root}/sylph"

    output:
    path("gtdb-r226-c200-dbv1.syldb"),                    emit: gtdb
    path("fungi-refseq-2025-10-11-c200-dbv1.syldb"),      emit: fungi
    path("imgvr_c200_v0.3.0.syldb"),                      emit: viral

    when:
    !file("${params.db_root}/sylph/gtdb-r226-c200-dbv1.syldb").exists() ||
    !file("${params.db_root}/sylph/fungi-refseq-2025-10-11-c200-dbv1.syldb").exists() ||
    !file("${params.db_root}/sylph/imgvr_c200_v0.3.0.syldb").exists()

    script:
    // Sylph databases are hosted at CMU (faust.compbio.cs.cmu.edu) with a
    // Google Cloud Storage mirror. We try both before failing.
    """
    set -euo pipefail

    fetch() {
        local fname="\$1"
        if [ -s "\$fname" ]; then
            echo "[setup_db] \$fname already present, skipping"
            return 0
        fi
        echo "[setup_db] Downloading \$fname..."
        wget -c --tries=3 --timeout=120 \\
            "http://faust.compbio.cs.cmu.edu/sylph-stuff/\$fname" \\
            || wget -c --tries=3 --timeout=120 \\
                "https://storage.googleapis.com/sylph-stuff/\$fname"
    }

    fetch gtdb-r226-c200-dbv1.syldb
    fetch fungi-refseq-2025-10-11-c200-dbv1.syldb
    fetch imgvr_c200_v0.3.0.syldb
    """
}

process DOWNLOAD_RESFINDER {
    storeDir "${params.db_root}/kma_resfinder/resfinder_db"

    output:
    path("all.fsa"), emit: done

    when:
    !file("${params.db_root}/kma_resfinder/resfinder_db/all.fsa").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Cloning ResFinder DB (~3 MB)..."
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git tmp_db
    cat tmp_db/*.fsa > all.fsa
    cp tmp_db/phenotypes.txt . 2>/dev/null || true
    cp tmp_db/config        . 2>/dev/null || true
    rm -rf tmp_db
    """
}

process SETUP_FASTQSCREEN {
    storeDir "${params.db_root}/fastqscreen"

    output:
    path("fastq_screen.conf"), emit: conf
    path("genomes/"),          emit: genomes

    when:
    !file("${params.db_root}/fastqscreen/fastq_screen.conf").exists()

    script:
    def fqs_url = params.fastqscreen_genomes_url ?: ''
    """
    set -euo pipefail
    mkdir -p genomes

    if [ -n "${fqs_url}" ]; then
        echo "[setup_db] Downloading FastQ Screen genome panel..."
        wget -c --tries=5 --timeout=120 "${fqs_url}" -O fastqscreen_genomes.tar.gz
        tar xzf fastqscreen_genomes.tar.gz -C genomes/ --strip-components=0
        rm -f fastqscreen_genomes.tar.gz
    else
        echo "[setup_db] WARNING: params.fastqscreen_genomes_url is not set."
        echo "[setup_db] FastQ Screen panel must be configured manually."
        echo "[setup_db] Place reference .fasta files under: ${params.db_root}/fastqscreen/genomes/"
    fi

    # Generate .fa symlinks and gzipped variants required by minimap2
    cd genomes
    for f in *.fasta; do
        [ ! -f "\$f" ] && continue
        base="\${f%.fasta}"
        ln -sf "\$f" "\${base}.fa" 2>/dev/null || true
        if [ ! -f "\${base}.fa.gz" ]; then
            gzip -c "\$f" > "\${base}.fa.gz"
        fi
    done
    cd ..

    # Auto-generate fastq_screen.conf from the available genomes
    {
        echo "# FastQ Screen config — generated automatically by EpiTaxMAG"
        echo "# Metagenomic surveillance panel"
        echo "ALIGNER   minimap2"
        echo "THREADS   8"
        echo ""
        for f in genomes/*.fasta; do
            [ ! -f "\$f" ] && continue
            name=\$(basename "\$f" .fasta)
            # Absolute prefix path; FastQ Screen appends .fa/.fa.gz automatically
            printf "DATABASE\\t%s\\t%s/genomes/%s\\n" "\$name" "\$(pwd)" "\$name"
        done
    } > fastq_screen.conf

    n=\$(ls genomes/*.fasta 2>/dev/null | wc -l)
    echo "[setup_db] FastQ Screen config generated with \$n genomes"
    """
}


// ──────────────────────────────────────────────────────────────────────
// MAG-phase databases
// ──────────────────────────────────────────────────────────────────────

process DOWNLOAD_CHECKM2 {
    storeDir "${params.db_root}/checkm2"

    output:
    path("CheckM2_database/uniref100.KO.1.dmnd"), emit: done

    when:
    params.run_assembly &&
    !file("${params.db_root}/checkm2/CheckM2_database/uniref100.KO.1.dmnd").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading CheckM2 DIAMOND DB (~1.7 GB)..."
    wget -c --tries=5 --timeout=120 \\
        https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
    tar xzf checkm2_database.tar.gz
    rm -f checkm2_database.tar.gz
    """
}

process DOWNLOAD_GTDBTK {
    storeDir "${params.db_root}/gtdbtk_r232"

    output:
    path("release232"), emit: done

    when:
    params.run_assembly && !params.skip_gtdbtk &&
    !file("${params.db_root}/gtdbtk_r232/release232").exists()

    script:
    // Pinned r232 download for reproducibility. The Australian mirror only
    // hosts the version-pinned tarball under the European GTDB host.
    """
    set -euo pipefail
    echo "[setup_db] Downloading GTDB-Tk r232 reference data (~60 GB compressed)..."
    wget -c --tries=5 --timeout=300 \\
        https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
    tar xzf gtdbtk_r232_data.tar.gz
    rm -f gtdbtk_r232_data.tar.gz
    test -d release232
    """
}

process DOWNLOAD_GENOMAD {
    storeDir "${params.db_root}/genomad"

    output:
    path("genomad_db"), emit: done

    when:
    params.run_assembly && !file("${params.db_root}/genomad/genomad_db").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading geNomad DB v1.9 (~840 MB)..."
    wget -c --tries=5 --timeout=120 \\
        https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz
    tar xzf genomad_db_v1.9.tar.gz
    rm -f genomad_db_v1.9.tar.gz
    test -d genomad_db
    """
}

process SETUP_AMRFINDERPLUS {
    // Uses the conda env containing amrfinder so the tool can fetch its DB.
    conda "${projectDir}/envs/nf-mag-amrfinderplus.yml"
    storeDir "${params.db_root}/amrfinderplus"

    output:
    path("latest/"), emit: db

    when:
    params.run_assembly && !file("${params.db_root}/amrfinderplus/latest").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading AMRFinderPlus DB..."
    mkdir -p latest
    amrfinder_update --database latest --force_update
    """
}

process SETUP_BAKTA {
    // Uses the conda env containing bakta so its tooling can fetch the DB.
    conda "${projectDir}/envs/nf-mag-bakta.yml"
    storeDir "${params.db_root}/bakta"

    output:
    path("db/"), emit: db

    when:
    params.run_assembly && !file("${params.db_root}/bakta/db").exists()

    script:
    def db_type = params.bakta_db_type ?: 'full'
    """
    set -euo pipefail
    echo "[setup_db] Downloading Bakta DB (${db_type}, ~84 GB for full / ~4 GB for light)..."
    bakta_db download --output db --type ${db_type}
    """
}

process DOWNLOAD_SOURMASH_GTDB {
    storeDir "${params.db_root}/sourmash"

    output:
    path("gtdb-rs226-reps-k31.zip"),  emit: db
    path("gtdb-rs226.lineages.csv"),  emit: lineages

    when:
    params.run_assembly && params.skip_gtdbtk &&
    !file("${params.db_root}/sourmash/gtdb-rs226-reps-k31.zip").exists()

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading Sourmash GTDB RS226 reps (k=31, ~3.6 GB)..."
    wget -c --tries=5 --timeout=120 \\
        https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip \\
        -O gtdb-rs226-reps-k31.zip
    echo "[setup_db] Downloading Sourmash GTDB RS226 lineages (~99 MB)..."
    wget -c --tries=5 --timeout=120 \\
        https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.csv
    """
}


// ──────────────────────────────────────────────────────────────────────
// Setup workflow — invoked from main.nf unless --skip_db_setup
// ──────────────────────────────────────────────────────────────────────

workflow SETUP_DB {

    main:

    // TAX databases (always required when --run_tax)
    DOWNLOAD_KRAKEN2()
    DOWNLOAD_KAIJU()
    DOWNLOAD_SYLPH()
    DOWNLOAD_RESFINDER()
    SETUP_FASTQSCREEN()

    // MAG databases (only when --run_assembly)
    if (params.run_assembly) {
        DOWNLOAD_CHECKM2()
        DOWNLOAD_GENOMAD()
        SETUP_AMRFINDERPLUS()
        SETUP_BAKTA()

        if (!params.skip_gtdbtk) {
            DOWNLOAD_GTDBTK()
        } else {
            DOWNLOAD_SOURMASH_GTDB()
        }
    }
}
