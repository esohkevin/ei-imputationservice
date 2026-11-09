#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
include {
  getVCF;
  getVcfs;
  getShapeitGeneticMap;
} from "${projectDir}/modules/phasing_and_imputation.nf"

workflow {
  getVCF().view()
  getVcfs().view()
  getShapeitGeneticMap().view()
}
*/


workflow {
  println "\n        IMPUTATION TEST"
  println ""
  println "idat_dir        = ${params.idat_dir}"
  println "manifest_bpm    = ${params.manifest_bpm}"
  println "manifest_csv    = ${params.manifest_csv}"
  println "cluster_file    = ${params.cluster_file}"

  if(params.build_ver == 'hg19') {
    println "fasta_ref      = ${params.fasta_ref}"
  }
  else {
    println "fasta_ref      = ${params.fasta_ref}"
    println "bam_alignment  = ${params.bam_alignment}"
  }

  println "output_prefix   = ${params.output_prefix}"
  println "output_dir      = ${params.output_dir}"
  println "containers_dir  = PATH WHERE CONTAINERS ARE STORED"
  println "account         = ${params.account}        # CHANGE TO YOURS"
  println "partition       = ${params.queue}      # CHANGE TO YOURS"
  println ""
  
  plink();
	
}

workflow.onComplete { 
  println "Workflow completed at: ${workflow.complete}"
  println "     Execution status: ${ workflow.success ? 'OK' : 'failed'}"
}

workflow.onError{
  println "workflow execution stopped with the following message: ${workflow.errorMessage}"
}

process plink() {
  tag "processing ... ${params.idat_dir}"
  label 'plink2'
  label 'idat_to_gtc'
  publishDir path: "${params.output_dir}/output"
  debug true
  
  script:
    """
    plink2 \
      --help \
      --file
    """
}

