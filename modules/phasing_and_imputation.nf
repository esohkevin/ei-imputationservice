def getVcfs() {
    return channel
               .fromFilePairs( params.vcf_dir + '*.vcf.{gz,gz.tbi}' )
               .ifEmpty { error: println "\nERROR: No VCF file found!\n" }
                  .map { vcfname, vcfindex -> tuple(vcfindex.first(), vcfindex.last()) }
}

def getVCF() {
    return channel.fromPath( params.vcf )
                  .ifEmpty { error: println "\nERROR: No VCF file found!\n" }
                  //.map { vcfname, vcfindex -> tuple(vcfindex.first(), vcfindex.last()) }
}

def listChromosomes() {
    return channel.of(1..22, 'X')
}

def getChromosomes() {
    if(params.build_ver == 'hg38') {
        if(params.autosome == true) {
           channel.of(1..22)
               .map { chr -> "chr${chr}" }
        } else {
           channel.of(1..22, 'X', 'X_PAR1', 'X_nonPAR', 'X_PAR2')
               .map { chr -> "chr${chr}" }
        }
    } else {
        if(params.autosome == true) {
           channel.of(1..22)
        } else {
           //channel.of(1..22, 'X', 'X_PAR1', 'X_nonPAR', 'X_PAR2')
           channel.of(1..22, 'X')
        }
    }
}

//-------------------------------------------------- REFERENCE PANELS ---------------------------------------------------------------

def getPanelMsav() {
    if(params.build_ver == 'hg38') {
    return channel.fromPath( params.panel_dir + "/" + params.panel + "/msav/*.msav" )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the MSAV reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { msav ->
                         tuple( msav.simpleName, msav)
                  }
    } else {
    return channel.fromPath( params.panel_dir + "/" + params.panel + "/msav/*.msav" )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the MSAV reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { msav ->
                         tuple( msav.simpleName.replaceAll(/chr/,"").toUpperCase(), msav)
                  }
    }
}

def getPanelBcf() {
    if(params.build_ver == 'hg38') {
    return channel
               .fromFilePairs( params.panel_dir + params.panel + "/bcf/chr*.{bcf,bcf.csi}", size: 2 )
               .ifEmpty { error: println "\nAn error occurred! Please check that the BCF reference files exist in ${params.panel_dir}/${params.panel}/\n" }
               .map { chr, ref_fileset ->
                      tuple(ref_fileset.first().simpleName, ref_fileset[0], ref_fileset[1])
               }
    } else {
    return channel
               .fromFilePairs( params.panel_dir + params.panel + "/bcf/chr*.{bcf,bcf.csi}", size: 2 )
               .ifEmpty { error: println "\nAn error occurred! Please check that the BCF reference files exist in ${params.panel_dir}/${params.panel}/\n" }
               .map { chr, ref_fileset ->
                      tuple(ref_fileset.first().simpleName.replaceFirst(/chr/,"").toUpperCase(), ref_fileset[0], ref_fileset[1])
               }
    }
}

def getPanelVcf() {
    if(params.build_ver == 'hg38') {
    return channel.fromFilePairs( params.panel_dir + params.panel + "/vcf/chr*.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the VCF reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { chr, ref_fileset ->
                      tuple(ref_fileset.first().simpleName, ref_fileset[0], ref_fileset[1])
                  }
    } else {
    return channel.fromFilePairs( params.panel_dir + params.panel + "/vcf/chr*.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the VCF reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { chr, ref_fileset ->
                      tuple(ref_fileset.first().simpleName.replaceFirst(/chr/,"").toUpperCase(), ref_fileset[0], ref_fileset[1])
                  }
    }
}

def getPanelBref() {
    if(params.build_ver == 'hg38') {
    return channel.fromPath( params.panel_dir + params.panel + "/bref/chr*.bref3" )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the BREF3 reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { ref ->
                      tuple(ref.simpleName, ref)
                  }
    } else {
    return channel.fromPath( params.panel_dir + params.panel + "/bref/chr*.bref3" ) 
                  .ifEmpty { error: println "\nAn error occurred! Please check that the BREF3 reference files exist in ${params.panel_dir}/${params.panel}/\n" }
                  .map { ref ->
                      tuple(ref.simpleName.replaceFirst(/chr/,"").toUpperCase(), ref)
                  }
    }
}


//----------------------------------------------------------------------------------------------------------------------------------

def getEagleHapmapGeneticMap() {
    return channel.fromPath( params.map_dir + "/eagle/genetic_map_hg19_withX.txt.gz" )
}

def getShapeitGeneticMap() {
    if(params.build_ver == 'hg38') {
    return channel
            .fromPath( params.map_dir + "/shapeit/" + "chr*.b38.gmap.gz" )
            .ifEmpty { error: println "\nAn error occurred! Please check that the map files exist in ${params.map_dir}/shapeit/\n" }
            .map { map -> tuple(map.simpleName, map) }
    } else {
    return channel
            .fromPath( params.map_dir + "/shapeit/" + "chr*.b37.gmap.gz" )
            .ifEmpty { error: println "\nAn error occurred! Please check that the map files exist in ${params.map_dir}/shapeit/\n" }
            .map { gmap -> tuple(gmap.simpleName.replaceAll(/chr/,"").toUpperCase(), gmap) }
    }
}

def getHapmapGeneticMap() {
    return channel
            .fromFilePairs( params.map_dir + "/hapmap/genetic_map_chr*_combined.txt", size: 1 )
            .ifEmpty { error: println "\nAn error occurred! Please check that the map files exist in ${params.map_dir}/hapmap/\n" }
            .map { chr, map_file ->
                   tuple( chr.replaceFirst(/genetic_map_chr/,"").toUpperCase(), map_file.first() )
            }
}

def getPlinkGeneticMap() {
    return channel
            .fromFilePairs( params.map_dir + '/plink/plink.chr*.GRCh37.map', size: 1)
            .ifEmpty { error: println "\nAn error occurred! Please check that the map files exist in ${params.map_dir}/plink/\n" }
            .map { chr, map_file ->
                   tuple( chr.replaceFirst(/^plink\.chr/,"").toUpperCase(), map_file.first() )
            }
}

def getPhasingRefPanel() {
    if(params.phase_tool.toUpperCase() == "BEAGLE5") {
        phaseref = getPanelBref()
    } else {
        phaseref = getPanelBcf()
    }
    return phaseref
}

def getImputationRefPanel() {
    if(params.impute_tool.toUpperCase() == 'MINIMAC4' ) {
        imputeref = getPanelMsav()
    } else if(params.impute_tool.toUpperCase() == "BEAGLE5") {
        imputeref = getPanelBref()
    } else {
        imputeref = getPanelBcf()
    }
    return imputeref
}


/*
process getPalindromicSnvs() {
    tag 'FILTER PALINDROMES'
    label 'Rbase'
    label 'extractPalindromes'
    input:
        tuple val(plink_base), path(plink_fileset)
    output:
        //publishDir path: "${params.outputDir}"
        path "${params.cohortName}.palindromic.snvs.txt"
    script:
        template "filterPalindromicVariants.r"
}
*/

process getContigsFromVcf() { 
    tag "processing ${vcf}..."
    label 'bcftools'
    label 'bigMemory'
    input:
        tuple \
            path(vcf), \
            path(index)
    output:
        tuple \
            path("contigs.list.txt"), \
            path(vcf), \
            path(index)
    script:
        """
        bcftools \
           reheader \
           -f ${params.fastaRef}.fai \
           ${vcf} | \
           tee tmp.vcf.gz | \
        bcftools \
           index \
           -ft \
           --threads ${task.cpus} \
           -o tmp.vcf.gz.tbi

        bcftools \
            index \
            --stats \
            tmp.vcf.gz | \
        awk '{print \$1,\$2}' > contigs.list.txt
        """
}

process getVcfGenomicIntervals() {
    tag "extracting invervals from VCF file..."
    label 'bcftools'
    cache 'lenient'
    input:
        tuple \
            path(contigs), \
            path(vcf), \
            path(index)
    output:
        path "*.bed"
    script:
        template 'get_gvcfinterval_list.sh'
}

process getChromFromVcf() {
    tag "processing ${vcf}..."
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            path(contigs), \
            path(vcf), \
            path(index)
    output:
        tuple \
            path("*.txt"), \
            path(vcf), \
            path(index)
    script:
        """
        for chrom in \$(awk '{print \$1}' ${contigs}); do 
            echo \${chrom} \
            > \${chrom}.txt; 
        done
        """
}

process indexVcf() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'longRun'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf}"), \
            path("${input_vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process getChromFromIntervalList() {
    tag "processing ${interval}..."
    input:
        path(interval)
    output:
        tuple \
            path("*.chrom"), \
            path("${interval}")
    script:
        """
        echo ${interval} > \$(awk '{print \$1}' ${interval}).chrom
        """
}

process splitVcfPerInterval() {
    tag "processing ${interval.baseName}..."
    label 'bcftools'
    label 'bigMemory'
    input:
        tuple \
            val(chrom), \
            path(vcf), \
            path(index), \
            path(interval)
    output:
        path("${interval.baseName}.vcf.gz*")
    script:
        """
        bcftools \
            view \
            -r \$(cat ${interval} | sed 's/ /:/1; s/ /-/1') \
            --threads ${task.cpus} \
            -Oz \
            -o ${interval.baseName}.vcf.gz \
            ${vcf}

        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${interval.baseName}.vcf.gz
        """
}

process validateVcfChunks() {
    tag "Validating VCF chunks..."
    input:
        path(vcf)
    output:
        path("*.validated.vcf.gz*")
    script:
        """
        readlink *.gz > vcflist.txt
        rm -rf *.vcf.gz*
        for vcf in \$(cat vcflist.txt); do
            vcfcontent=\$(zgrep -v '^#' \${vcf} | head -1 | cut -f1)
            if [[ ! \${vcfcontent} == "" ]]; then
                ln -s \${vcf} \$(basename \${vcf/.vcf.gz/.validated.vcf.gz})
                ln -s \${vcf}.tbi \$(basename \${vcf/.vcf.gz/.validated.vcf.gz}.tbi)
            fi
        done
        """
}

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf}"), \
            path("${input_vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process removeDupVarsAndIndexVcf() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf.simpleName}.nodups.vcf.gz"), \
            path("${input_vcf.simpleName}.nodups.vcf.gz.tbi")
    script:
        """
        #bcftools \
        #    sort \
        #    -Oz \
        #    -o ${input_vcf.baseName}-sorted.vcf.gz \
        #    ${input_vcf}

        bcftools \
            norm \
            -m +any \
            -D \
            --threads ${task.cpus} \
            -Oz \
            ${input_vcf} | \
        bcftools \
            view \
            --threads ${task.cpus} \
            -m2 \
            -M2 \
            -Oz | \
        tee ${input_vcf.simpleName}.nodups.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${input_vcf.simpleName}.nodups.vcf.gz.tbi \
        """
}

process splitVcfByChrom() {
    tag "processing chr${chrom}"
    label 'plink'
    label 'splitVCF'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}.vcf.gz")
    script:
        """
        plink2 \
            --vcf ${input_vcf} \
            --chr ${chrom} \
            --vcf-half-call missing \
            --export vcf-4.2 bgz id-paste='iid' \
            --threads ${task.cpus} \
            --out "chr${chrom}"
        """
}

process sortIndexVcf() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf.baseName}-sorted.vcf.gz"), \
            path("${input_vcf.baseName}-sorted.vcf.gz.tbi")
    script:
        """
        bcftools \
            sort \
            --threads ${task.cpus} \
            -Oz \
            -o ${input_vcf.baseName}-sorted.vcf.gz \
            ${input_vcf}

        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf.baseName}-sorted.vcf.gz
        """
}

process checkStrand() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'alignGenotypes'
    publishDir \
        pattern: "chr${chrom}-aligned.vcf.gz", \
        path: "${params.output_dir}/aligned/vcf/", \
        mode: 'copy'
    publishDir \
        pattern: "chr${chrom}-aligned.log", \
        path: "${params.output_dir}/aligned/logs/", \
        mode: 'copy'
    input:
	tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}-aligned.vcf.gz"), \
            path("chr${chrom}-aligned.log")
    script:
    if(params.exclude_sample == "NULL")
        """
        conform-gt \
            gt="${input_vcf}" \
            ref="${ref_vcf}" \
            chrom=${chrom} \
            match=POS \
            out="chr${chrom}-aligned"
        """
    else
	"""
	conform-gt \
	    gt="${input_vcf}" \
	    ref="${ref_vcf}" \
	    chrom=${chrom} \
	    match=POS \
            excludesamples="${params.exclude_sample}" \
	    out="chr${chrom}-aligned"
	"""
}

process concatVcfs() {
    tag "concatenating VCFs into ${params.output_prefix}_aligned.vcf.gz..."
    label 'bcftools'
    label 'mediumMemory'
    publishDir \
        path: "${params.output_dir}/aligned/combined/", \
        mode: 'copy'
    input:
        path(vcf)
    output:
        path("${params.output_prefix}_aligned.vcf.gz*")
    script:
        """
        ls *.vcf.gz | sort -V > aligned_vcf_list.txt
        bcftools \
            concat \
            -f aligned_vcf_list.txt \
            -a \
            --threads ${task.cpus} \
            -Oz | \
        bcftools \
            norm \
            -m+ \
            --threads ${task.cpus} \
            -Oz | \
        bcftools \
            view \
            -v snps \
            -m2 -M2 \
            --threads ${task.cpus} \
            -Oz | \
        tee ${params.output_prefix}_aligned.vcf.gz | \
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            -o ${params.output_prefix}_aligned.vcf.gz.tbi
        """
}

process beaglephase() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'phaseGenotypes'
    input:
	tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
	publishDir path: "${params.output_dir}", mode: 'copy'
	tuple \
            val(chrom), \
            path("chr${chrom}.*.vcf.gz")
    script:
    if(params.impute == true)
	"""
	beagle \
	    gt="${input_vcf}" \
	    ref="${ref_vcf}" \
	    map="${geneticMap}" \
	    chrom=${chrom} \
	    burnin="${params.burnit}" \
	    iterations="${params.mainit}" \
	    ne=15000 \
	    impute=${params.impute} \
	    nthreads=${task.cpus} \
	    out="chr${chrom}.dose
	"""
    else
        """
        beagle \
            gt=${input_vcf} \
            ref=${ref_vcf} \
            map=${geneticMap} \
            chrom=${chrom} \
            burnin=${params.burnit} \
            iterations=${params.mainit} \
            ne=15000 \
            impute=${params.impute} \
            nthreads=${task.cpus} \
            out="chr${chrom}.phased
        """
}

process beaglePhaseWithoutRef() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'phaseGenotypes'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}.*.vcf.tmp.gz")
    script:
        """
        beagle \
            gt=${input_vcf} \
            chrom=${chrom} \
            burnin=${params.burnit} \
            iterations=${params.mainit} \
            ne=15000 \
            nthreads=${task.cpus} \
            out=chr${chrom}.phased.tmp
        """
}

process eaglePhaseWithoutRef() {
    tag "processing chr${chrom}"
    label 'eagle'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(geneticMap)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}.${params.output_prefix}.eagle2.noref.tmp.vcf.gz"), \
            path("chr${chrom}.${params.output_prefix}.eagle2.noref.log")            
    script:
        """
        eagle \
          --vcf=${input_vcf} \
          --geneticMapFile=${geneticMap} \
          --chrom=${chrom} \
          --numThreads=${task.cpus} \
          --Kpbwt=${params.kpbwt} \
          --vcfOutFormat=z \
          --outPrefix=chr${chrom}.${params.output_prefix}.eagle2.tmp.noref 2>&1 | \
          tee chr${chrom}.${params.output_prefix}.eagle2.noref.log
        """
}

process eaglePhaseWithRef() {
    tag "processing chr${chrom}"
    label 'eagle'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}.${params.output_prefix}.eagle2.tmp.vcf.gz"), \
            path("chr${chrom}.${params.output_prefix}.eagle2.log")
    script:
        """
        eagle \
          --vcfTarget=${input_vcf} \
          --vcfRef=${ref_vcf} \
          --geneticMapFile=${geneticMap} \
          --chrom=${chrom} \
          --numThreads=${task.cpus} \
          --Kpbwt=${params.kpbwt} \
          --vcfOutFormat=z \
          --outPrefix=chr${chrom}.${params.output_prefix}.eagle2.tmp 2>&1 | \
          tee chr${chrom}.${params.output_prefix}.eagle2.log
        """
}

process shapeitPhaseWithoutRef() {
    tag "processing chromosome ${chrom}"
    label 'shapeit'
    label 'phase_shapeit'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(geneticMap)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf.baseName}.shapeit4.noref.tmp.vcf.gz")
    script:
    if(params.wgs == true)
        """
        shapeit4 \
          --input ${input_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --mcmc-iterations 10b,1p,2b,1p,1b,1p,1b,1p,10m \
          --pbwt-depth ${params.pbwt} \
          --sequencing \
          --log ${input_vcf.baseName}.shapeit4.noref.log \
          --output ${input_vcf.baseName}.shapeit4.noref.tmp.vcf.gz
        """
    else
        """
        shapeit4 \
          --input ${input_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --mcmc-iterations 10b,1p,2b,1p,1b,1p,1b,1p,10m \
          --pbwt-depth ${params.pbwt} \
          --log ${input_vcf.baseName}.shapeit4.noref.log \
          --output ${input_vcf.baseName}.shapeit4.noref.tmp.vcf.gz
        """
}

process shapeitPhaseWithRef() {
    tag "processing chr${chrom}"
    label 'shapeit'
    label 'phase_shapeit'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
        tuple \
            val(chrom), \
            path("${input_vcf.baseName}.shapeit4.vcf.gz")
    script:
    if(params.wgs == true)
        """
        shapeit4 \
          --input ${input_vcf} \
          --reference ${ref_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --sequencing \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --log ${input_vcf.baseName}.shapeit4.log \
          --output ${input_vcf.baseName}.shapeit4.vcf.gz
        """
    else
        """
        shapeit4 \
          --input ${input_vcf} \
          --reference ${ref_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --log ${input_vcf.baseName}.shapeit4.log \
          --output ${input_vcf.baseName}.shapeit4.vcf.gz
        """
}

process updateVcfHeader() {
    tag "updating header for ${vcf}..."
    label 'bcftools'
    label 'bigMemory'
    publishDir \
        path: "${params.output_dir}/phased/", \
        mode: 'copy'
    input:
        tuple \
            val(chrom), \
            path(vcf)
    output:
        tuple \
            val(chrom), \
            path("${vcf.baseName}.fixed.vcf.gz")
    script:
        """
        bcftools \
            reheader \
            --fai ${params.fastaRef}.fai \
            -o ${vcf.baseName}.fixed.vcf.gz \
            ${vcf}
        """
}

process getm3vcf() {
    tag "processing chr${chrom}"
    label 'minimac3'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index) 
    output:
        publishDir path: "${params.panel_dir}/m3vcfs/", mode: 'copy'
        tuple \
            val(chrom), \
            path("chr${chrom}*")
    script:
        """
        Minimac3 \
          --refHaps ${input_vcf} \
          --processReference \
          --prefix chr${chrom} \
          --cpus ${task.cpus} \
        """
}

process createLegendFile() {
    tag "processing chr${chrom}"
    label 'vcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index)
    output:
        publishDir path: "${params.output_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path(input_vcf), \
            path("${input_vcf.simpleName}.legend.gz")
    script:
        """
        echo "id position a0 a1 all.aaf" > header
        vcftools \
            --gzvcf ${input_vcf} \
            --freq \
            --out ${chrom}
        sed 's/:/\\t/g' ${chrom}.frq | \
            sed 1d | \
            awk '{print \$1":"\$2" "\$2" "\$5" "\$7" "\$8}' \
            > ${input_vcf.simpleName}.legend
        cat header ${input_vcf.simpleName}.legend | \
            bgzip \
            > ${input_vcf.simpleName}.legend.gz
        """
}

process prepareChrXPanel() {
    tag "processing chr${chrom}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        publishDir path: "${params.output_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${input_vcf.simpleName}_PAR1.vcf.gz"), \
            path("${input_vcf.simpleName}_PAR1,nonPAR.vcf.gz.tbi"), \
            path("${input_vcf.simpleName}_PAR2.vcf.gz"), \
            path("${input_vcf.simpleName}_PAR2,nonPAR.vcf.gz.tbi"), \
            path("${input_vcf.simpleName}_nonPAR.vcf.gz"), \
            path("${input_vcf.simpleName}_nonPAR.vcf.gz.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads \
            ${task.cpus} \
            ${input_vcf}
        bcftools \
            view \
            -r X:60001-2699520 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_PAR1.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_PAR1.vcf.gz.tbi
        bcftools \
            view \
            -r X:154931044-155260560 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_PAR2.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_PAR2.vcf.gz.tbi
        bcftools \
            view \
            -r X:2699521-154931043 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_nonPAR.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_nonPAR.vcf.gz.tbi            
        """
}

process imputeVariantsWithMinimac4Old() {
    tag "processing chr${chrom}"
    label 'minimac4_old'
    label 'mem_minimac4'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(erate), \
            path(ref_vcf), \
            path(rec), \
            path(target_vcf), \
            path(target_vcf_index)
    output:
        publishDir path: "${params.output_dir}/${params.panel}/", mode: 'copy'
        tuple \
            val(chrom), \
            path("chr${chrom}.dose.vcf.gz"), \
            path("chr${chrom}.info.gz"), \
            path("chr${chrom}.logfile")
    script:
        """
        minimac4 \
          --refHaps ${ref_vcf} \
          --haps ${target_vcf} \
          --log chr${chrom} \
          --prefix chr${chrom} \
          --cpus ${task.cpus} \
          --allTypedSites \
          --ChunkLengthMb 5 \
          --ChunkOverlapMb 0.0025 \
          --noPhoneHome

        bgzip -f chr${chrom}.info
        """
}

process imputeVariantsWithMinimac4() {
    tag "processing chr${chrom}"
    label 'minimac4'
    label 'mem_minimac4'
    //cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(msav), \
            path(target_vcf), \
            path(target_vcf_index)
    output:
        publishDir path: "${params.output_dir}/${params.panel}/", mode: 'copy'
        tuple \
            val(chrom), \
            path("*${chrom}.dose.vcf.gz"), \
            path("*${chrom}.empiricalDose.vcf.gz")
    script:
    if(params.build_ver == "hg19")
        """
        minimac4 \
          --all-typed-sites \
          --chunk ${params.chunk_size} \
          --threads ${task.cpus} \
          --format GT,HDS,DS,GP \
          --output chr${chrom}.dose.vcf.gz \
          --empirical-output chr${chrom}.empiricalDose.vcf.gz \
          ${msav} \
          ${target_vcf}
        """
    else
        """
        minimac4 \
          --all-typed-sites \
          --chunk ${params.chunk_size} \
          --threads ${task.cpus} \
          --format GT,HDS,DS,GP \
          --output ${chrom}.dose.vcf.gz \
          --empirical-output ${chrom}.empiricalDose.vcf.gz \
          ${msav} \
          ${target_vcf}
        """

}

process collectIntervalsPerChromosome() {
    tag "Collecting intervals per chromosome..."
    label 'bcftools'
    label 'gatkVariantCaller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    //storeDir "${params.output_dir}/vcf/"
    input:
        path(vcfs)
    output:
        path("*_vcfs_list.txt")
    script:
        """
        ls *.{vcf.gz,bcf} | awk '{print \$1,\$1}' > vcfs_list.txt

        while read line; do 
            data=( \$line );
            chr=\$(basename \${data[0]} | sed 's/_/ /g; s/\\./ /1' | awk '{print \$1,\$2}')
            echo \$chr \$(readlink \${data[1]})
        done < vcfs_list.txt | sort -V -k1 -k2 > vcf_chr_list.txt


        for chrom in \$(awk '{print \$1}' vcf_chr_list.txt | sort -V | uniq); do
            awk -v chr="\${chrom}" '\$1 == chr' vcf_chr_list.txt | sort -V -k1 -k2 > \${chrom}_vcfs_list.txt
        done
        """
}

process concatPerChromIntervalVcfs() {
    tag "Concatenating VCF files per chromosome..."
    label 'bcftools'
    label 'gatkVariantCaller'
    //storeDir "${params.output_dir}/vcf/"
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'copy'
    input:
        path(vcf_list)
    output:
        tuple \
            path("*.tmp.txt"), \
            path("*_${params.output_prefix}.vcf.gz"), \
            path("*_${params.output_prefix}.vcf.gz.tbi")
    script:
        """
        chrom=\$(awk '{print \$1}' ${vcf_list} | uniq)
        awk '{print \$3}' ${vcf_list} > \${chrom}_concat.list

        for vcf in \$(cat \${chrom}_concat.list); do
            [[ ! -e "\${vcf}.csi"  ]] && \
            bcftools \
                index \
                --threads ${task.cpus} \
                \${vcf}
        done

        bcftools \
            concat \
            -a \
            --threads ${task.cpus} \
            -Oz \
            -f \${chrom}_concat.list | \
        tee \${chrom}_${params.output_prefix}.vcf.gz | \
        bcftools \
            index \
            --threads ${task.cpus} \
            -ft \
            -o \${chrom}_${params.output_prefix}.vcf.gz.tbi

        echo \$chrom > \${chrom}.tmp.txt
        """
}

process concatPerChromosomeVcfs() {
    tag "Concatenating all VCF files into ${params.output_prefix}.vcf.gz..."
    label 'bcftools'
    label 'gatkVariantCaller'
    publishDir \
        path: "${params.output_dir}/vcf/", \
        mode: 'move'
    input:
        tuple \
            path("*.tmp.txt"), \
            path("*_${params.output_prefix}.vcf.gz"), \
            path("*_${params.output_prefix}.vcf.gz.tbi")
    output:
        path("${params.output_prefix}.vcf.{gz,gz.tbi}")
    script:
        """
        readlink *.gz | sort -V > concat.list

        bcftools \
            concat \
            -a \
            --threads ${task.cpus} \
            -Oz \
            -f concat.list | \
        tee ${params.output_prefix}.vcf.gz | \
        bcftools \
            index \
            --threads ${task.cpus} \
            -ft \
            -o ${params.output_prefix}.vcf.gz.tbi
        """
}

