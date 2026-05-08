#!/usr/bin/env nextflow

/*
 * SETUP_DB — Automatic download of pipeline databases.
 *
 * Each process uses `storeDir` so files persist across runs and across
 * pipeline restarts. Nextflow skips a process when all its outputs are
 * already in the storeDir, so the workflow is idempotent without an
 * explicit `when:` guard. Use `--skip_db_setup` to bypass entirely.
 *
 * The workflow emits a single `ready` channel that downstream subworkflows
 * (TAX, MAG) gate on, so no taxonomy/AMR process starts before its
 * databases finish downloading.
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

    // Declare every file Kraken2/Bracken expect so storeDir keeps them all.
    // Without `database*.kmer_distrib`, Bracken would fail at runtime.
    output:
    path("hash.k2d"),                         emit: hash
    path("opts.k2d"),                         emit: opts
    path("taxo.k2d"),                         emit: taxo
    path("inspect.txt"),                      emit: inspect, optional: true
    path("ktaxonomy.tsv"),                    emit: tax,     optional: true
    path("seqid2taxid.map"),                  emit: seqid,   optional: true
    path("database*.kmer_distrib"),           emit: bracken, optional: true

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

    // The 2024-08-14 tarball ships the .fmi alongside the NCBI taxonomy
    // dump. Declare all three so storeDir keeps them — without nodes.dmp
    // and names.dmp Kaiju and kaiju2table fail at runtime.
    output:
    path("kaiju_db_refseq_ref.fmi"), emit: fmi
    path("nodes.dmp"),               emit: nodes
    path("names.dmp"),               emit: names

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading Kaiju RefSeq ref index 2024-08-14 (~36 GB compressed)..."
    wget -c --tries=5 --timeout=120 \\
        https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz
    tar xzf kaiju_db_refseq_ref_2024-08-14.tgz
    rm -f kaiju_db_refseq_ref_2024-08-14.tgz
    """
}

process DOWNLOAD_SYLPH {
    storeDir "${params.db_root}/sylph"

    output:
    path("gtdb-r226-c200-dbv1.syldb"),                    emit: gtdb
    path("fungi-refseq-2025-10-11-c200-dbv1.syldb"),      emit: fungi
    path("imgvr_c200_v0.3.0.syldb"),                      emit: viral

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

    script:
    def fqs_url   = params.fastqscreen_genomes_url ?: ''
    def final_dir = "${params.db_root}/fastqscreen/genomes"
    """
    set -e
    shopt -s nullglob
    mkdir -p genomes

    if [ -n "${fqs_url}" ]; then
        echo "[setup_db] Downloading FastQ Screen genome panel..."
        wget -c --tries=5 --timeout=120 "${fqs_url}" -O fastqscreen_genomes.tar.gz
        # The Zenodo tarball nests fastas under a top-level `genomes/` directory.
        # Extract into a temp dir, then flatten everything into our genomes/.
        mkdir -p _extract
        tar xzf fastqscreen_genomes.tar.gz -C _extract/
        find _extract -name '*.fasta' -exec mv {} genomes/ \\;
        rm -rf _extract fastqscreen_genomes.tar.gz
    else
        echo "[setup_db] WARNING: params.fastqscreen_genomes_url is not set."
        echo "[setup_db] FastQ Screen panel must be configured manually."
        echo "[setup_db] Place reference .fasta files under: ${final_dir}/"
    fi

    # Generate .fa symlinks and gzipped variants required by minimap2
    pushd genomes >/dev/null
    for f in *.fasta; do
        base="\${f%.fasta}"
        ln -sf "\$f" "\${base}.fa" 2>/dev/null || true
        if [ ! -f "\${base}.fa.gz" ]; then
            gzip -c "\$f" > "\${base}.fa.gz"
        fi
    done
    popd >/dev/null

    # Auto-generate fastq_screen.conf with absolute paths to the FINAL storeDir
    # location, so it works once Nextflow moves the outputs.
    {
        echo "# FastQ Screen config — generated automatically by EpiTaxMAG"
        echo "# Metagenomic surveillance panel"
        echo "ALIGNER   minimap2"
        echo "THREADS   8"
        echo ""
        pushd genomes >/dev/null
        for f in *.fasta; do
            name="\${f%.fasta}"
            # FastQ Screen appends .fa/.fa.gz automatically
            printf "DATABASE\\t%s\\t%s/%s\\n" "\$name" "${final_dir}" "\$name"
        done
        popd >/dev/null
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

process DOWNLOAD_SKANI {
    // skani v0.3 ships a single-file sketch DB for GTDB r226. Bluenote-1577
    // hosts the pre-built archive on the CMU compbio server. The folder
    // extracted by tar (`skani_gtdb_r226-v0.3/`) is what `skani search -d`
    // expects as its input.
    storeDir "${params.db_root}/skani"

    output:
    path("skani_gtdb_r226-v0.3"), emit: db

    script:
    """
    set -euo pipefail
    echo "[setup_db] Downloading skani GTDB r226 sketch DB (~38 GB compressed)..."
    wget -c --tries=5 --timeout=300 \\
        http://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz
    tar xzf skani_gtdb_r226-v0.3.tar.gz
    rm -f skani_gtdb_r226-v0.3.tar.gz
    test -d skani_gtdb_r226-v0.3
    """
}


// ──────────────────────────────────────────────────────────────────────
// Setup workflow — invoked from main.nf unless --skip_db_setup.
// Emits a single `ready` channel that downstream subworkflows gate on.
// ──────────────────────────────────────────────────────────────────────

workflow SETUP_DB {

    main:

    // TAX databases (always required when --run_tax)
    DOWNLOAD_KRAKEN2()
    DOWNLOAD_KAIJU()
    DOWNLOAD_SYLPH()
    DOWNLOAD_RESFINDER()
    SETUP_FASTQSCREEN()

    // Collect TAX-phase done signals
    ch_tax_done = DOWNLOAD_KRAKEN2.out.hash
        .mix(
            DOWNLOAD_KAIJU.out.fmi,
            DOWNLOAD_SYLPH.out.gtdb,
            DOWNLOAD_RESFINDER.out.done,
            SETUP_FASTQSCREEN.out.conf
        )
        .collect()

    if (params.run_assembly) {
        DOWNLOAD_CHECKM2()
        DOWNLOAD_GENOMAD()
        SETUP_AMRFINDERPLUS()
        SETUP_BAKTA()

        // Download only the DB for the selected MAG taxonomy tool.
        switch (params.taxonomy_tool) {
            case 'skani':
                DOWNLOAD_SKANI()
                ch_mag_tax = DOWNLOAD_SKANI.out.db
                break
            case 'sourmash':
                DOWNLOAD_SOURMASH_GTDB()
                ch_mag_tax = DOWNLOAD_SOURMASH_GTDB.out.db
                break
            case 'gtdbtk':
                DOWNLOAD_GTDBTK()
                ch_mag_tax = DOWNLOAD_GTDBTK.out.done
                break
            default:
                error "Unknown taxonomy_tool: '${params.taxonomy_tool}'. Use 'skani', 'sourmash' or 'gtdbtk'."
        }

        ch_mag_done = DOWNLOAD_CHECKM2.out.done
            .mix(
                DOWNLOAD_GENOMAD.out.done,
                SETUP_AMRFINDERPLUS.out.db,
                SETUP_BAKTA.out.db,
                ch_mag_tax
            )
            .collect()

        ch_ready = ch_tax_done.mix(ch_mag_done).collect().map { 'ready' }
    } else {
        ch_ready = ch_tax_done.map { 'ready' }
    }

    emit:
    ready = ch_ready
}
