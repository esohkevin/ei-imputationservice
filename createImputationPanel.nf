#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getVcf;
    splitVcfByChrom;
    beaglephase;
    eaglePhaseWithoutRef;
    getEagleHapmapGeneticMap;
    getVcfIndex;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
    getValidationExitStatus;
    createLegendFile;
    getm3vcf;
    prepareChrXPanel;
} from "${projectDir}/modules/custom_panel.nf"

include {
    jobCompletionMessage;
} from "${projectDir}/includes/inspection.nf"


workflow {
    chromosome = getChromosomes()
    vcf = getPhasedVcf()
    validateVcf(vcf).map { chr, vcf_file -> tuple(chr.baseName, vcf_file) }.set { chrom_vcf }
    getVcfIndex(chrom_vcf).view().set { vcf_fileset }

    createLegendFile(vcf_fileset).view()
    m3vcf = getm3vcf(vcf_fileset).view()
    if(params.autosome == false) {
        channel.of('X')
            .join()
        prepareChrXPanel(vcf_fileset).view()
    }

}

workflow.onComplete {
    msg = jobCompletionMessage()
    if(params.email == 'NULL') {
        println msg
    } 
    else {
        println msg
        sendMail(
            //from: 'eshkev001@myuct.ac.za',
            to: params.email,
            subject: '[genemapis status]',
            body: msg
        )
    }
}
