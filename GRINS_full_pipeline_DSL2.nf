#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Combined GRINS detection workflow, modernized for Nextflow DSL2.
 *
 * Replaces:
 *   1. detect_in_genome.nf
 *   2. manual asmash_genbanks symlink step
 *   3. GRINS_detection_from_BOWTIE.nf
 *
 * Main improvements:
 *   - DSL2 workflow/process calling syntax.
 *   - path/val tuple inputs instead of DSL1 `from` input clauses.
 *   - named process outputs with `emit`.
 *   - container/profile support through nextflow.config labels.
 *   - stageInMode 'copy' and cp -L for container-safe file staging.
 *   - optional reuse of already generated antiSMASH output with --antismash_indir.
 */

params.indir = 'genomes_fasta'
params.outdir = 'full_run'
params.sensitivity = 'sensitive'
params.w_size = 150
params.s_size = 30
params.min_size = 500
params.antismash_cpus = 8
params.bowtie2_cpus = 1
params.with_plots = 'yes'
params.antismash_bin = '/Users/kian/grins/bin/antismash'

/*
 * Optional directory containing already completed antiSMASH outputs.
 *
 * Leave blank to run antiSMASH normally.
 *
 * Use this to skip antiSMASH:
 *
 *   --antismash_indir full_run/antismash
 *
 * For an input FASTA named:
 *
 *   genomes_fasta/sequence.fasta
 *
 * this expects either:
 *
 *   full_run/antismash/sequence/sequence.json
 *   full_run/antismash/sequence/sequence.gbk
 *
 * or:
 *
 *   full_run/antismash/sequence.json
 *   full_run/antismash/sequence.gbk
 */
params.antismash_indir = ''

process RUN_ANTISMASH {
    label 'run_antismash'
    cpus params.antismash_cpus
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/antismash",
        mode: 'copy',
        saveAs: { filename -> filename == 'output' ? acc : null }

    input:
    tuple val(acc), path(genome_file)

    output:
    tuple val(acc), path('output'), emit: antismash_dir

    script:
    """
    set -euo pipefail

    GENOME_INPUT="${genome_file}"
    GENOME_FOR_ANTISMASH="${acc}.fasta"

    if [ "\$GENOME_INPUT" != "\$GENOME_FOR_ANTISMASH" ]; then
        cp -f "\$GENOME_INPUT" "\$GENOME_FOR_ANTISMASH"
    fi

    ${params.antismash_bin} \
        -c ${task.cpus} \
        --taxon bacteria \
        --output-dir output \
        --verbose \
        --genefinding-tool prodigal \
        "\$GENOME_FOR_ANTISMASH"
    """
}

process STAGE_EXISTING_ANTISMASH {
    label 'py3'
    tag "$acc"
    stageInMode 'copy'

    input:
    tuple val(acc), path(json_file), path(gbk_file)

    output:
    tuple val(acc), path('output'), emit: antismash_dir

    script:
    """
    set -euo pipefail

    mkdir -p output

    cp -L "${json_file}" "output/${acc}.json"
    cp -L "${gbk_file}" "output/${acc}.gbk"

    if [[ ! -s "output/${acc}.json" ]]; then
        echo "ERROR: Existing antiSMASH JSON file is missing or empty after staging: output/${acc}.json" >&2
        exit 1
    fi

    if [[ ! -s "output/${acc}.gbk" ]]; then
        echo "ERROR: Existing antiSMASH GenBank file is missing or empty after staging: output/${acc}.gbk" >&2
        exit 1
    fi
    """
}

process ANTISMASH2GFF3 {
    label 'py3'
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/antismash.gff3",
        mode: 'copy'

    input:
    tuple val(acc), path(antismash_dir)

    output:
    tuple val(acc), path("${acc}.gff3"), emit: regions_gff3

    script:
    """
    set -euo pipefail

    python3 "${workflow.projectDir}/antismash2gff3.py" \
        --input "${antismash_dir}/${acc}.json" \
        --output "${acc}.gff3"
    """
}

process COLLECT_ANTISMASH_GBK {
    label 'py3'
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/asmash_genbanks",
        mode: 'copy'

    input:
    tuple val(acc), path(antismash_dir)

    output:
    tuple val(acc), path("${acc}.gbk"), emit: gbk

    script:
    """
    set -euo pipefail

    if [[ ! -f "${antismash_dir}/${acc}.gbk" ]]; then
        echo "ERROR: Expected antiSMASH GenBank file not found: ${antismash_dir}/${acc}.gbk" >&2
        echo "Available files in antiSMASH output:" >&2
        find "${antismash_dir}" -maxdepth 2 -type f -print >&2
        exit 1
    fi

    cp -L "${antismash_dir}/${acc}.gbk" "${acc}.gbk"
    """
}

process SPLIT_IN_WINDOWS {
    label 'py3'
    tag "$acc"
    stageInMode 'copy'

    input:
    tuple val(acc), path(genome_file)

    output:
    tuple val(acc), path("${acc}.fasta"), path("${acc}_windows.fasta"), emit: windows

    script:
    """
    set -euo pipefail

    echo "split_seq_into_windows pipe-header-v2" >&2

    python3 "${workflow.projectDir}/split_seq_into_windows.py" \
        --input "${acc}.fasta" \
        --format fasta \
        --w_size ${params.w_size} \
        --s_size ${params.s_size} \
        --output "${acc}_windows.fasta"
    """
}

process BOWTIE2_ALIGN {
    label 'bowtie2'
    cpus params.bowtie2_cpus
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/bam",
        mode: 'copy'

    input:
    tuple val(acc), path(fasta), path(windows)

    output:
    tuple val(acc), path("${acc}.bam"), path(fasta), emit: bam

    script:
    def sensFlag = params.sensitivity == 'very-sensitive' ? '--very-sensitive' : '--sensitive'
    """
    set -euo pipefail

    mkdir -p idx
    bowtie2-build "${fasta}" "idx/${acc}"

    bowtie2 \
        -f \
        --end-to-end \
        ${sensFlag} \
        -a \
        --time \
        --threads ${task.cpus} \
        -x "idx/${acc}" \
        -U "${windows}" | \
    samtools view -b -o "${acc}.bam" -
    """
}

process MERGE_BAM_WINDOWS {
    label 'py3'
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/duplicated.gff3",
        mode: 'copy'

    input:
    tuple val(acc), path(bam), path(fasta)

    output:
    tuple val(acc), path("${acc}.duplicated.gff3"), emit: dup_gff3

    script:
    """
    set -euo pipefail

    python3 "${workflow.projectDir}/produce_windows_from_bam.py" \
        --input "${bam}" \
        --output "${acc}.duplicated.gff3" \
        --w_size ${params.w_size} \
        --min_size ${params.min_size}
    """
}

process INTERSECT_ASMASH_DUPS {
    label 'bedtools'
    tag "$acc"
    stageInMode 'copy'

    publishDir "${params.outdir}/bgcdups.gff3",
        mode: 'copy'

    input:
    tuple val(acc), path(dup_gff3), path(regions_gff3)

    output:
    tuple val(acc), path("${acc}.bgcdups.gff3"), emit: bgcdups

    script:
    """
    set -euo pipefail

    if [ ! -s "${dup_gff3}" ]; then
        printf "##gff-version 3\\n" > "${acc}.bgcdups.gff3"
        exit 0
    fi

    bedtools intersect -wo -a "${dup_gff3}" -b "${regions_gff3}" | \\
    awk -F '\\t' 'BEGIN { OFS="\\t" } {
        gsub(/ID=region/, "Region=region", \$9);
        \$9 = \$9 ";" \$18;
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9;
    }' > "${acc}.bgcdups.gff3"

    if [ ! -s "${acc}.bgcdups.gff3" ]; then
        printf "##gff-version 3\\n" > "${acc}.bgcdups.gff3"
    fi
    """
}

process GRINSPRED {
    label 'py3'
    tag "$genome"
    stageInMode 'copy'

    publishDir "${params.outdir}/grins/outputs",
        mode: 'copy',
        pattern: 'output',
        saveAs: { filename -> filename == 'output' ? genome : null }

    publishDir "${params.outdir}/grins/res",
        mode: 'copy',
        pattern: 'GRINS_detected_in_genomes_and_BGCs.txt',
        saveAs: { filename -> "${genome}.txt" }

    input:
    tuple val(genome), path(gbk_file), path(dup_gff3)

    output:
    tuple val(genome), path('output'), path('GRINS_detected_in_genomes_and_BGCs.txt'), emit: grins

    script:
    """
    set -euo pipefail

    mkdir -p gbk dup output/genomes_GRINS output/GRINS.gff3 output/GRINS_BGC.gff3 output/plots

    cp -L "${gbk_file}" "gbk/${genome}.gbk"
    cp -L "${dup_gff3}" "dup/${genome}.duplicated.gff3"

    python3 "${workflow.projectDir}/GRINS_detection_from_BOWTIE.py" \
        --seq_input gbk \
        --dupl_input dup \
        --seq_output output/genomes_GRINS \
        --GRINS_output output/GRINS.gff3 \
        --GRINS_BGC_output output/GRINS_BGC.gff3 \
        --plot_output output/plots \
        --with_plots ${params.with_plots}
    """
}

workflow {
    if (!(params.sensitivity in ['sensitive', 'very-sensitive'])) {
        error "Invalid --sensitivity '${params.sensitivity}'. Use 'sensitive' or 'very-sensitive'."
    }

    genomes_ch = Channel
        .fromPath("${params.indir}/*", type: 'file', checkIfExists: true)
        .filter { it.name ==~ /.+\.(fa|fasta|fna|fas)$/ }
        .ifEmpty { error "No FASTA files found in ${params.indir}. Expected .fa, .fasta, .fna, or .fas files." }
        .map { genomefa ->
            def acc = genomefa.name.replaceFirst(/\.(fa|fasta|fna|fas)$/, '')
            tuple(acc, genomefa)
        }

    SPLIT_IN_WINDOWS(genomes_ch)

    def antismash_dir_ch

    if (params.antismash_indir != null && params.antismash_indir.toString().trim() != '') {
        log.info "Using existing antiSMASH results from: ${params.antismash_indir}"
        log.info "Skipping RUN_ANTISMASH process."

        existing_antismash_files_ch = genomes_ch.map { acc, genomefa ->
            def root = file(params.antismash_indir)

            if (!root.exists()) {
                error "--antismash_indir does not exist: ${params.antismash_indir}"
            }

            def nested_json = file("${params.antismash_indir}/${acc}/${acc}.json")
            def nested_gbk = file("${params.antismash_indir}/${acc}/${acc}.gbk")

            def flat_json = file("${params.antismash_indir}/${acc}.json")
            def flat_gbk = file("${params.antismash_indir}/${acc}.gbk")

            if (nested_json.exists() && nested_gbk.exists()) {
                tuple(acc, nested_json, nested_gbk)
            }
            else if (flat_json.exists() && flat_gbk.exists()) {
                tuple(acc, flat_json, flat_gbk)
            }
            else {
                error """
Could not find existing antiSMASH files for genome '${acc}'.

Expected one of these layouts:
  ${params.antismash_indir}/${acc}/${acc}.json
  ${params.antismash_indir}/${acc}/${acc}.gbk

or:
  ${params.antismash_indir}/${acc}.json
  ${params.antismash_indir}/${acc}.gbk

Make sure the antiSMASH file basename matches the FASTA basename in --indir.
"""
            }
        }

        STAGE_EXISTING_ANTISMASH(existing_antismash_files_ch)
        antismash_dir_ch = STAGE_EXISTING_ANTISMASH.out.antismash_dir
    }
    else {
        RUN_ANTISMASH(genomes_ch)
        antismash_dir_ch = RUN_ANTISMASH.out.antismash_dir
    }

    ANTISMASH2GFF3(antismash_dir_ch)
    COLLECT_ANTISMASH_GBK(antismash_dir_ch)

    BOWTIE2_ALIGN(SPLIT_IN_WINDOWS.out.windows)
    MERGE_BAM_WINDOWS(BOWTIE2_ALIGN.out.bam)

    INTERSECT_ASMASH_DUPS(
        MERGE_BAM_WINDOWS.out.dup_gff3.join(ANTISMASH2GFF3.out.regions_gff3)
    )

    GRINSPRED(
        COLLECT_ANTISMASH_GBK.out.gbk.join(MERGE_BAM_WINDOWS.out.dup_gff3)
    )
}
