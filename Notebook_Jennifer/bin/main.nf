#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMsg() {
  log.info """
  USAGE:
  The typical command for running the pipeline is as follows:
  nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" -profile nova
  nextflow run main.nf --genome GENOME.fasta --reads_file READS_PATH.txt -profile nova
  
  Mandatory arguments:
    --genome             Genome fasta file, against which reads will be mapped
    --reads              Paired end reads in fastq.gz format

  Optional configuration arguments:
    --gmap_build_app     Link to gmap_build executable
    --gsnap_app          Link to gsnap executable
    --samtools_app       Link to samtools executable
    --featureCounts      Link to featureCounts executable
  """
}

/* ==== Define Processes ==== */
process index_reference {
  tag "$fasta_gz"
  label "gmap_build"
  publishDir "${params.outdir}/ref"

  input:   // genome.fasta.gz
  path fasta_gz

  output: // genome_name, gmapdb directory
  tuple val("${fasta_gz.simpleName}"), path("${params.outdir}/ref")

  script:
  """
  #! /usr/bin/env bash
  gmap_build -d "${fasta_gz.simpleName}" -D "${params.outdir}/ref" --gunzip ${fasta_gz}
  """
}

/* ==== Main Workflow === */
workflow {
  channel.fromPath(params.genome) | view
  channel.fromPath(params.reads) | view

  genome_ch = channel.fromPath(params.genome) | index_reference
}

