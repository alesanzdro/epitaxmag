#!/usr/bin/env nextflow

/*
 * SETUP_DB — Descarga automatica de bases de datos
 * Solo descarga lo que falta. Si todo existe, no hace nada.
 */

process DOWNLOAD_KRAKEN2 {
    storeDir "${params.db_root}/k2_pluspf_20251015"

    output:
    path("hash.k2d"), emit: done

    when:
    !file("${params.db_root}/k2_pluspf_20251015/hash.k2d").exists()

    script:
    """
    echo "Descargando Kraken2 PlusPF (~100 GB)..."
    wget -c -q https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz
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
    echo "Descargando Kaiju RefSeq (~34 GB)..."
    wget -c -q https://bioinformatics-centre.github.io/kaiju/downloads/kaiju_db_refseq_ref.tgz
    tar xzf kaiju_db_refseq_ref.tgz
    rm -f kaiju_db_refseq_ref.tgz
    """
}

process DOWNLOAD_SYLPH {
    storeDir "${params.db_root}/sylph"

    output:
    path("gtdb-r226-c200-dbv1.syldb"), emit: done

    when:
    !file("${params.db_root}/sylph/gtdb-r226-c200-dbv1.syldb").exists()

    script:
    """
    echo "Descargando Sylph GTDB r226 (~14 GB)..."
    wget -c -q https://storage.googleapis.com/sylph-stuff/v0.7-databases/gtdb-r226-c200-dbv1.syldb
    echo "Descargando Sylph Fungi (~700 MB)..."
    wget -c -q https://storage.googleapis.com/sylph-stuff/v0.7-databases/fungi-refseq-2024-07-25-c200-v0.3.syldb
    echo "Descargando Sylph Viral (~2 GB)..."
    wget -c -q https://storage.googleapis.com/sylph-stuff/v0.7-databases/imgvr_c200_v0.3.0.syldb
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
    echo "Descargando ResFinder DB..."
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git tmp_db 2>/dev/null || true
    cat tmp_db/*.fsa > all.fsa
    cp tmp_db/phenotypes.txt . 2>/dev/null || true
    rm -rf tmp_db
    """
}

process DOWNLOAD_CHECKM2 {
    storeDir "${params.db_root}/checkm2"

    output:
    path("CheckM2_database/uniref100.KO.1.dmnd"), emit: done

    when:
    params.run_assembly && !file("${params.db_root}/checkm2/CheckM2_database/uniref100.KO.1.dmnd").exists()

    script:
    """
    echo "Descargando CheckM2 (~3 GB)..."
    wget -c -q https://zenodo.org/records/5571251/files/checkm2_database.tar.gz
    tar xzf checkm2_database.tar.gz
    rm -f checkm2_database.tar.gz
    """
}

process DOWNLOAD_GTDBTK {
    storeDir "${params.db_root}/gtdbtk_r232"

    output:
    path("release232"), emit: done

    when:
    params.run_assembly && !params.skip_gtdbtk && !file("${params.db_root}/gtdbtk_r232/release232").exists()

    script:
    """
    echo "Descargando GTDB-Tk r232 (~100 GB comprimido)..."
    wget -c -q https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    tar xzf gtdbtk_data.tar.gz
    rm -f gtdbtk_data.tar.gz
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
    echo "Descargando geNomad (~1.4 GB)..."
    wget -c -q https://zenodo.org/records/8339387/files/genomad_db_v1.7.tar.gz
    tar xzf genomad_db_v1.7.tar.gz
    rm -f genomad_db_v1.7.tar.gz
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
    mkdir -p genomes

    if [ -n "${fqs_url}" ]; then
        echo "Descargando panel FastQ Screen..."
        wget -c -q "${fqs_url}" -O fastqscreen_genomes.tar.gz
        tar xzf fastqscreen_genomes.tar.gz -C genomes/ --strip-components=0
        rm -f fastqscreen_genomes.tar.gz
    else
        echo "AVISO: No se ha definido params.fastqscreen_genomes_url"
        echo "El panel de FastQ Screen debe configurarse manualmente."
        echo "Crear genomas .fasta en: ${params.db_root}/fastqscreen/genomes/"
    fi

    # Generar .fa symlinks y .fa.gz para FastQ Screen
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

    # Generar fastq_screen.conf automaticamente
    echo "# FastQ Screen config — generado automaticamente por EpiTaxMAG" > fastq_screen.conf
    echo "# Panel de vigilancia metagenomica" >> fastq_screen.conf
    echo "ALIGNER   minimap2" >> fastq_screen.conf
    echo "THREADS   8" >> fastq_screen.conf
    echo "" >> fastq_screen.conf

    for f in genomes/*.fasta; do
        [ ! -f "\$f" ] && continue
        name=\$(basename "\$f" .fasta)
        # Ruta absoluta al prefijo (sin extension, FastQ Screen añade .fa/.fa.gz)
        echo "DATABASE\t\${name}\t\$(pwd)/genomes/\${name}" >> fastq_screen.conf
    done

    echo "FastQ Screen config generado con \$(ls genomes/*.fasta 2>/dev/null | wc -l) genomas"
    """
}

process SETUP_AMRFINDERPLUS {
    storeDir "${params.db_root}/amrfinderplus"

    output:
    path("latest/"), emit: db

    when:
    params.run_assembly && !file("${params.db_root}/amrfinderplus/latest").exists()

    script:
    """
    echo "Descargando AMRFinderPlus DB..."
    amrfinder_update --database latest 2>/dev/null || \
    amrfinder --update --database latest 2>/dev/null || \
    echo "AVISO: No se pudo descargar AMRFinderPlus DB automaticamente"
    """
}

process SETUP_BAKTA {
    storeDir "${params.db_root}/bakta"

    output:
    path("db/"), emit: db

    when:
    params.run_assembly && !file("${params.db_root}/bakta/db").exists()

    script:
    def db_type = params.bakta_db_type ?: 'full'
    """
    echo "Descargando Bakta DB ${db_type}..."
    bakta_db download --output db --type ${db_type} 2>/dev/null || \
    echo "AVISO: No se pudo descargar Bakta DB automaticamente."
    echo "Descargar manualmente: bakta_db download --output \${PWD}/db --type ${db_type}"
    """
}

process DOWNLOAD_SOURMASH_GTDB {
    storeDir "${params.db_root}/sourmash"

    output:
    path("gtdb-rs226-k31.zip"), emit: done

    when:
    params.run_assembly && params.skip_gtdbtk && !file("${params.db_root}/sourmash/gtdb-rs226-k31.zip").exists()

    script:
    """
    echo "Descargando Sourmash GTDB RS226 reps k31 (~3.7 GB)..."
    wget -c -q https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip -O gtdb-rs226-reps-k31.zip
    echo "Descargando lineages (~99 MB)..."
    wget -c -q https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.csv
    """
}

workflow SETUP_DB {
    main:

    // TAX databases
    DOWNLOAD_KRAKEN2()
    DOWNLOAD_KAIJU()
    DOWNLOAD_SYLPH()
    DOWNLOAD_RESFINDER()
    SETUP_FASTQSCREEN()

    // MAG databases
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
