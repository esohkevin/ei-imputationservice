#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
