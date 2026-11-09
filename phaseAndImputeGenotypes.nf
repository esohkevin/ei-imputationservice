#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVCF;
    getVcfs;
    getChromosomes;
    //listChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getEagleHapmapGeneticMap;
    getShapeitGeneticMap;
    getPanelMsav;
    getPanelBcf;
    getPanelVcf;
    getPanelBref;
    getPhasingRefPanel;
    getImputationRefPanel;
    //
    //getRefPanel;
    //getKgpPanel;
    //getGmv2p1;
    //getCustomVcfPanel;
    //getCustomM3Panel;
    //
    splitVcfByChrom;
    getVcfGenomicIntervals;
    splitVcfPerInterval;
    getChromFromIntervalList;
    validateVcfChunks;
    getChromFromVcf;
    getContigsFromVcf;
    collectIntervalsPerChromosome;
    concatPerChromIntervalVcfs;
    concatPerChromosomeVcfs;
    beaglephase;
    beaglePhaseWithoutRef;
    eaglePhaseWithoutRef;
    eaglePhaseWithRef;
    shapeitPhaseWithRef;
    shapeitPhaseWithoutRef;
    updateVcfHeader;
    getVcfIndex;
    getVcfIndex as getPhasedVcfIndex;
    removeDupVarsAndIndexVcf;
    imputeVariantsWithMinimac4;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
} from "${projectDir}/modules/custom_panel.nf"


workflow {

    getChromosomes()
        .set { chromosome }

    //getVCF()

    vcf = getVcfs()

    getContigsFromVcf(vcf)
        .set{ contigs }

    getVcfGenomicIntervals(contigs)
        .flatten()
        .map { interval ->
            tuple( "${interval.simpleName}", interval )
        }
        .groupTuple()
        .set{ chrom_interval }

    getChromFromVcf(contigs)
        .flatMap { chroms, vcf, index ->
            chroms.collect { chrom ->
                tuple("${chrom.simpleName}", vcf, index)
            }
        }
        .join(chrom_interval)
        .flatMap { chrom, vcf, index, intervals ->
            intervals.collect { interval ->
                tuple(chrom, vcf, index, interval)
            }
        }
        .set { chrom_vcf_interval }

    splitVcfPerInterval(chrom_vcf_interval)
        .collect()
        .set { vcfs }

    
    /****************************************************************************
    * First group tuple by 0 (chrom) to collect all intervals per chrom together
    * in order to merge with per chrom genetic maps. Then you can flatten the 
    * the merged set and then group again by 1 (key)
    *****************************************************************************/
    
    validateVcfChunks(vcfs)
        .flatten()
        .map { vcf ->
            tuple("${vcf.simpleName}", vcf )
        }
        .groupTuple(by:0)
        .set { vcfFileset }

    if(params.phase == true && params.impute == false) {
       println "\nMODE: PHASE ONLY\n"
       if(params.phase_tool == 'eagle2') {
            println "\nPHASING WITH EAGLE2\n"
            //if(params.with_ref == true) {
            //    geneticmap = getEagleHapmapGeneticMap()
            //    if(params.panel == "kgp") {
            //        refpanel = getKgpPanel() 
            //    } else {
            //        refpanel = getGmv2p1()
            //    }
            //    refpanel.combine(geneticmap).set { panel_map }
            //    vcfFileset.join(panel_map).set { phase_input }
            //    eaglePhaseWithRef(phase_input)
            //        .map { chrom, vcf, log -> tuple(chrom, vcf) }
            //        .set { ephased }
            //    updateVcfHeader(ephased)
            //        .set { phased }
            //} else {
            //    geneticmap = getEagleHapmapGeneticMap()
            //    vcfFileset.combine(geneticmap).set { phase_input }
            //    eaglePhaseWithoutRef(phase_input)
            //    eaglePhaseWithRef(phase_input)
            //        .map { chrom, vcf, log -> tuple(chrom, vcf) }
            //        .set { ephased }
            //    updateVcfHeader(ephased)
            //        .set { phased }
            //}
        } else if(params.phase_tool == 'shapeit4') {
            println "\nPHASING WITH SHAPEIT4\n"
            if(params.with_ref == true) {
                println "\nPHASING WITH REFERENCE\n"
                //geneticmap = getShapeitGeneticMap()
                ////refpanel = getKgpPanel()
                //refpanel = getGmv2p1()
                //refpanel.join(geneticmap).set { panel_map }
                //vcfFileset.join(panel_map).set { phase_input }
                //shapeitPhaseWithRef(phase_input)
                //    .set { sphased }
                //updateVcfHeader(sphased)
                //    .set { phased }
            } else {
                 geneticmap = getShapeitGeneticMap()
                 vcfFileset
                    .join(geneticmap)
                    .flatMap { chrom, interval_filesets, gmap ->
                        interval_filesets.collect { interval_fileset ->
                            tuple(chrom, "${interval_fileset.baseName.replaceAll(/.vcf*/,"").replaceAll(/.gz/,"")}", interval_fileset, gmap)
                        }
                    }
                    .groupTuple(by:1)
                    .map { chrom, key, vcffileset, gmap ->
                        tuple(chrom[0], vcffileset[0], vcffileset[1], gmap[0] )
                    }
                    .set { phase_input }
                shapeitPhaseWithoutRef(phase_input)
                    .set { sphased }
                updateVcfHeader(sphased)
                    .map { chr, vcf -> vcf }
                    .collect()
                    .set { phased_per_interval }
                collectIntervalsPerChromosome(phased_per_interval)
                    .flatten()
                    .set { vcfs_per_chrom_list }
                concatPerChromIntervalVcfs(vcfs_per_chrom_list)
                    .map { chr, vcf, index ->
                        tuple("${chr.simpleName}", vcf, index)
                    }
                    .set { phased_vcf_fileset }

                //geneticmap = getShapeitGeneticMap()
                //vcfFileset
                //    .join(geneticmap).view()
                //    .set { phase_input }
                //shapeitPhaseWithoutRef(phase_input)
                //    .set { sphased }
                //updateVcfHeader(sphased)
                //    .collect()
                //    .set { phased }
            }
        } else if(params.phase_tool == 'beagle5') {
            println "\nPHASING WITH BEAGLE5\n"
            if(params.with_ref == true) {
                println "\nPHASING WITH REFERENCE\n"
                //geneticmap = getPlinkGeneticMap()
                ////refpanel = getKgpPanel()
                //refpanel = getGmv2p1()
                //refpanel.join(geneticmap).set { panel_map }
                //vcfFileset.join(panel_map).set { phase_input }
                //beaglePhaseWithoutRef(phase_input)
                //    .set { bphased }
                //updateVcfHeader(bphased)
                //    .set { phased }
            } else {
                beaglePhaseWithoutRef(vcfFileset)
                    .set { bphased }
                updateVcfHeader(bphased)
                    .set { phased }
            }
        }

       //vcfs_per_chrom_list = collectIntervalsPerChromosome(phased).flatten()
       //vcfs_per_chrom = concatPerChromIntervalVcfs(vcfs_per_chrom_list).collect().view()
    }
    else if(params.phase == false && params.impute == true) {
        println "\nMODE: IMPUTE ONLY\n"
        
        //if(params.impute_tool == 'minimac4') {
        //    //vcf = getPhasedVcf()
        //    //validateVcf(vcf).map { chr, vcf_file, vcf_index -> tuple(chr.baseName, vcf_file, vcf_index) }.set { vcf_fileset }
        //    //getCustomM3Panel().set{ minimac_ref_panel }
        //    getRefPanel()
        //        .set{ minimac_ref_panel }
        //    minimac_ref_panel.join(vcfFileset).set{ minimac_input }
        //    imputeVariantsWithMinimac4(minimac_input)
        //} else { // Build commands for imputation using other tools and panels
        //    getRefPanel()
        //}
    } else if(params.phase == true && params.impute == true) {
        println "\nMODE: PHASE AND IMPUTE\n"

        if(params.phase_tool == 'eagle2' && params.impute_tool == 'minimac4') {
            if(params.with_ref == true) {
                getEagleHapmapGeneticMap().set { geneticmap }
                getPhasingRefPanel().set { refpanel }
                refpanel.combine(geneticmap).set { panel_map }
                vcfFileset.join(panel_map).set { phase_input }
                eaglePhaseWithRef(phase_input)
                    .map { chrom, vcf, log -> tuple(chrom, vcf) }
                    .set { ephased }
                updateVcfHeader(ephased)
                    .set { phased }
            } else {
                getEagleHapmapGeneticMap().set { geneticmap}
                vcfFileset.combine(geneticmap).set { phase_input }
                eaglePhaseWithoutRef(phase_input)
                    .map { chrom, vcf, log -> tuple(chrom, vcf) }
                    .set { ephased }
                updateVcfHeader(ephased)
                    .map { chr, vcf -> vcf }
                    .collect()
                    .set { phased_per_interval }
                collectIntervalsPerChromosome(phased_per_interval)
                    .flatten()
                    .set { vcfs_per_chrom_list }
                concatPerChromIntervalVcfs(vcfs_per_chrom_list)
                    .set { phased }
            }

            getPhasedVcfIndex(phased)
                .set { phased_vcf_fileset }
            getImputationRefPanel()
                .set{ minimac_ref_panel }
            minimac_ref_panel.join(phased_vcf_fileset)
                .set{ minimac_input }
            imputeVariantsWithMinimac4(minimac_input)

        } else if(params.phase_tool == 'shapeit4' && params.impute_tool == 'minimac4') {
            if(params.with_ref == true) {
                geneticmap = getShapeitGeneticMap()
                getPhasingRefPanel()
                    .join(geneticmap)
                    .set { panel_map }
                vcfFileset
                    .join(panel_map)
                    .flatMap { chrom, interval_filesets, ref, refindex, gmap ->
                        interval_filesets.collect { interval_fileset ->
                            tuple(chrom, "${interval_fileset.baseName.replaceAll(/.vcf*/,"").replaceAll(/.gz/,"")}", interval_fileset, ref, refindex, gmap)
                        }
                    }
                    .groupTuple(by:1)
                    .map { chrom, key, vcffileset, ref, refindex, gmap ->
                        tuple(chrom[0], vcffileset[0], vcffileset[1], ref[0], refindex[0], gmap[0] ) 
                    }
                    .set { phase_input }
                shapeitPhaseWithRef(phase_input)
                    .set { sphased }
                updateVcfHeader(sphased)
                    .map { chr, vcf -> vcf }
                    .collect()
                    .set { phased_per_interval }
                collectIntervalsPerChromosome(phased_per_interval)
                    .flatten()
                    .set { vcfs_per_chrom_list }
                concatPerChromIntervalVcfs(vcfs_per_chrom_list)
                    .map { chr, vcf, index -> 
                        tuple("${chr.simpleName}", vcf, index)
                    }
                    .set { phased_vcf_fileset }
            } else {
                geneticmap = getShapeitGeneticMap()
                vcfFileset
                    .join(panel_map)
                    .flatMap { chrom, interval_filesets, gmap ->
                        interval_filesets.collect { interval_fileset ->
                            tuple(chrom, "${interval_fileset.baseName.replaceAll(/.vcf*/,"").replaceAll(/.gz/,"")}", interval_fileset, gmap)
                        }
                    }
                    .groupTuple(by:1)
                    .map { chrom, key, vcffileset, gmap ->
                        tuple(chrom[0], vcffileset[0], vcffileset[1], gmap[0] )
                    }
                    .set { phase_input }
                shapeitPhaseWithRef(phase_input)
                    .set { sphased }
                updateVcfHeader(sphased)
                    .map { chr, vcf -> vcf }
                    .collect()
                    .set { phased_per_interval }
                collectIntervalsPerChromosome(phased_per_interval)
                    .flatten()
                    .set { vcfs_per_chrom_list }
                concatPerChromIntervalVcfs(vcfs_per_chrom_list)
                    .map { chr, vcf, index ->
                        tuple("${chr.simpleName}", vcf, index)
                    }
                    .set { phased_vcf_fileset }
            }

            getImputationRefPanel()
                .set{ minimac_ref_panel }
            minimac_ref_panel.join(phased_vcf_fileset)
                .set{ minimac_input }
            imputeVariantsWithMinimac4(minimac_input)

        } else if(params.phase_tool == 'beagle5') {
            refpanel = getPhasingRefPanel()
            geneticmap = getPlinkGeneticMap()
            vcfFileset
                .join(refpanel)
                .join(geneticmap)
                .set { beagle_input }
           beaglephase(beagle_input)
        }
    } else {
        error println "\nWORKFLOW STOPPED: Please select a run mode - 'phase' and/or 'impute' -\n"
    }
} 
