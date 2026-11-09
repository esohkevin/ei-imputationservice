#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getEagleHapmapGeneticMap;
    getShapeitGeneticMap;
    getPanelVcf;
    getVCF;
    splitVcfByChrom;
    sortIndexVcf;
    checkStrand;
    concatVcfs;
    getVcfIndex;
    removeDupVarsAndIndexVcf;
    getVcfIndex as getPhasedVcfIndex;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
} from "${projectDir}/modules/custom_panel.nf"

workflow {

    if(params.autosome == true) {
       getChromosomes()
          .filter(Number)
          .set { chromosome }
    } else {
       getChromosomes()
          .set { chromosome }
    }

    getVCF()
       .set { vcf }

    chromosome
       .combine(vcf)
       .set { split_vcf_input }

    getPanelVcf()
       .set { refpanel }

    //getBref()
    //    .set { refpanel }

    getHapmapGeneticMap()
       .set { geneticMap }

    splitVcfByChrom(split_vcf_input)
       .set { per_chr_vcf }

    removeDupVarsAndIndexVcf(per_chr_vcf)
       .set { vcf_fileset }

    vcf_fileset
       .map { chr, vcf, index -> tuple("${chr}", vcf, index) }
       .set { vcfFileset }

    vcfFileset
	.map { chrom, vcfFile, vcfIndex -> tuple( "${chrom}", vcfFile, vcfIndex ) }
	.join( refpanel )
	.join( geneticMap )
	.set { checkstrand_input }

    alignedVcfs = checkStrand( checkstrand_input ).map { chr, vcf, log -> tuple(chr, vcf) }
    aligned_indexed_vcfs = getVcfIndex(alignedVcfs).map { chr, vcf, index -> tuple(vcf, index) }.collect()
    concatVcfs(aligned_indexed_vcfs)
}
