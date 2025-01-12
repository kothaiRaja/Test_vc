nextflow.enable.dsl = 2

params.reads = "/home/kothai/cq-git-sample/test/data/test/STAR/trimmed_sample_*_Aligned.sortedByCoord.out.bam"
params.genome = "/home/kothai/cq-git-sample/test/data/test/genome.fa"
params.fasta_index = "/home/kothai/cq-git-sample/test/data/test/genome.fa.fai"
params.genome_dict = "/home/kothai/cq-git-sample/test/data/test/genome.dict"
params.outdir= "output"
params.filtered_vcf = "/home/kothai/cq-git-sample/test/data/test/merged.filtered.recode.vcf.gz"
params.filtered_vcf_index = "/home/kothai/cq-git-sample/test/data/test/merged.filtered.recode.vcf.gz.tbi"
params.denylist = "/home/kothai/cq-git-sample/test/data/test/denylist.bed"
params.scatter_count = 3 
params.mode = "test"
params.genomedb = 'GRCh38.86'
params.merge_vcf = false 
params.remove_duplicates = false
params.validation_stringency = 'STRICT'
params.vep_cache = "/home/kothai/cq-git-sample/test/data/test/vep_cache"
params.clinvar = "/home/kothai/cq-git-sample/test/data/test/clinvar/clinvar.vcf.gz"
params.clinvartbi = "/home/kothai/cq-git-sample/test/data/test/clinvar/clinvar.vcf.gz.tbi"



process SAMTOOLS_SORT_INDEX {
	tag { sample_id }
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/sorted_bam", mode: "copy"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam ${bam}

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam
    """
}

process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/filtered_bam", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai")

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -f 2 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam

    # Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
    """
}



process SAMTOOLS_FLAGSTAT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), path("${sample_id}_stats_report.txt")

    script:
    """
    # Validate BAM file to ensure it is not corrupted
    samtools quickcheck -v ${sorted_bam} || (echo "BAM file validation failed for ${sample_id}" && exit 1)

    # Run samtools flagstat to generate alignment metrics
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt

    # Generate additional quality metrics with samtools stats
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt
    """
}



process GATK_MARK_DUPLICATES {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/mark_duplicates", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bam_index)

    output:
    tuple val(sample_id), 
          path("${sample_id}_marked_duplicates.bam"), 
          path("${sample_id}_marked_duplicates.bai"), 
          path("${sample_id}_dup_metrics.txt")

    script:
    """
    gatk MarkDuplicates \
        -I ${sorted_bam} \
        -O ${sample_id}_marked_duplicates.bam \
        -M ${sample_id}_dup_metrics.txt \
        --CREATE_INDEX true \
		--REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \
        --VALIDATION_STRINGENCY ${params.validation_stringency ?: 'LENIENT'}

    """
}

process SPLIT_NCIGAR_READS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/split_ncigar", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(genome_fasta), path (index), path (genome_dict)

    output:
    tuple val(sample_id), 
          path("${sample_id}_split.bam"), 
          path("${sample_id}_split.bai")

    script:
    """
    gatk SplitNCigarReads \
        -R ${genome_fasta} \
        -I ${bam} \
        -O ${sample_id}_split.bam \
        --create-output-bam-index true  
		
	"""
}

process SAMTOOLS_CALMD {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), 
          path(bam), 
          path(bai) 
          path(genome_fasta) 
		  path (index)
		  
    output:
    tuple val(sample_id), 
          path("${sample_id}_calmd.bam"), 
          path("${sample_id}_calmd.bam.bai") 
         

    script:
    """
    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${sample_id}_calmd.bam

    # Index the updated BAM file
    samtools index ${sample_id}_calmd.bam

    """
}


process GATK_RECALIBRATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recalibrated_bams", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(genome_fasta), path(index), path(dict), path(known_variants), path(known_variants_index)

    output:
    tuple val(sample_id), 
          path("${sample_id}_recalibrated.bam"), 
          path("${sample_id}_recalibrated.bai"), 
          path("${sample_id}_recal_data.table")

    script:
    """
    # Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -R ${genome_fasta} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal_data.table 
		

    # Step 2: ApplyBQSR
    gatk ApplyBQSR \
        -R ${genome_fasta} \
        -I ${bam} \
        --bqsr-recal-file ${sample_id}_recal_data.table \
        -O ${sample_id}_recalibrated.bam

	# Step 3: Validation
	gatk ValidateSamFile \
		-I ${sample_id}_recalibrated.bam \
		-MODE SUMMARY 


    """
}

process BED_TO_INTERVAL_LIST {
    tag { bed_file.name }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/intervals", mode: "copy"

    input:
    path bed_file
    path genome_fasta
    path genome_dict

    output:
    path("${bed_file.baseName}.interval_list")

    script:
    """
    gatk BedToIntervalList \
        -I ${bed_file} \
        -O ${bed_file.baseName}.interval_list \
        -SD ${genome_dict}
    """
}

process SCATTER_INTERVAL_LIST {
    tag "Scatter interval list"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/scattered_intervals", mode: 'copy'

    input:
    path interval_list
    path genome_dict

    output:
    path "*.interval_list"

    script:
    """
    mkdir -p scattered_intervals
    gatk IntervalListTools \
        --INPUT ${interval_list} \
        --OUTPUT scattered_intervals \
        --SCATTER_COUNT ${params.scatter_count} \
        --UNIQUE true

    # Move and rename the output files to the working directory with unique names
    for f in scattered_intervals/*/*; do
        dir_name=\$(basename \$(dirname "\$f"))  # Get subdirectory name
        file_name=\$(basename "\$f")            # Get original file name
        mv "\$f" "scattered_intervals/\${dir_name}_\${file_name}.interval_list"
    done

    # Move the uniquely renamed files to the working directory
    mv scattered_intervals/*.interval_list .
    rm -r scattered_intervals
	
	# Log the generated intervals
    echo "Scattered intervals:"
    cat *.interval_list
    """
}





process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }

     container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi") 

    script:
	def intervals_args = interval_list.collect { "--intervals ${it}" }.join(' ')
	
    """
    gatk HaplotypeCaller \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${genome} \
    --output output_${sample_id}.vcf.gz \
    -I $bam \
    --standard-min-confidence-threshold-for-calling 5.0 \
    --dont-use-soft-clipped-bases true \
    --min-base-quality-score 10 \
	--output-mode EMIT_VARIANTS_ONLY \
    ${intervals_args} \
    --verbosity DEBUG

    """
}

process BCFTOOLS_STATS {
    tag "bcftools_stats"

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h8b25389_0" 
    publishDir "${params.outdir}/bcftools_stats", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)   

    output:
    tuple val(sample_id), path("haplotypecaller_stats_${sample_id}.txt")   

    script:
    """
    # Generate stats
    bcftools stats output_${sample_id}.vcf.gz > haplotypecaller_stats_${sample_id}.txt

    """
}


process GATK_VARIANT_FILTER {
    tag "${sample_id}_variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi")

    script:
    """
    gatk VariantFiltration \
        -R ${genome} \
        -V ${vcf_file} \
        --cluster-window-size 35 \
        --cluster-size 3 \
        --filter-name "LowQual" --filter-expression "QUAL < 20.0" \
        --filter-name "LowQD" --filter-expression "QD < 1.5" \
        --filter-name "HighFS" --filter-expression "FS > 60.0" \
        --filter-name "LowMQ" --filter-expression "MQ < 30.0" \
        --filter-name "HighSOR" --filter-expression "SOR > 4.0" \
        --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < -5.0" \
        --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < -3.0" \
        -O ${sample_id}_filtered.vcf.gz

    # Validate the filtered VCF
    if [ ! -s ${sample_id}_filtered.vcf.gz ]; then
        echo "Error: Filtered VCF is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}


process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
    publishDir "${params.outdir}/variant_summary", mode: "copy"

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index)

    output:
    tuple val(sample_id), path("filtered_variants_summary_${sample_id}.txt")

    script:
    """
    # Generate a summary of filtered variants
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${filtered_vcf} > filtered_variants_summary_${sample_id}.txt

    # Validate the output
    if [ ! -s filtered_variants_summary_${sample_id}.txt ]; then
        echo "Error: Summary of filtered variants is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}


process ANNOTATE_INDIVIDUAL_VARIANTS {
    tag "${sample_id}_annotate"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index)
	path (snpEffJar)
	path (snpEffConfig)
	path (snpEffDbDir)
	val (genomedb)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"), path("${sample_id}.annotated.summary.html")

    script:
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${filtered_vcf} > ${sample_id}.annotated.vcf
		



    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats ${sample_id}.annotated.summary.html \
        ${filtered_vcf} > /dev/null
    """
}

process ANNOTATE_INDIVIDUAL_VARIANTS_VEP {
    tag "${sample_id}_vep_annotate"

    container "https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index)  // Input VCF and its index
    path vep_cache                                                  // VEP cache directory
    path clinvar_vcf                                                // ClinVar VCF file
    path clinvar_index                                              // ClinVar Tabix index file

    output:
    tuple val(sample_id), path("${sample_id}.vep.annotated.vcf"), path("${sample_id}.vep.summary.html")

    script:
    """
    # Annotate using Ensembl VEP with ClinVar
    vep --input_file ${filtered_vcf} \
        --output_file ${sample_id}.vep.annotated.vcf \
        --stats_file ${sample_id}.vep.summary.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN

    # Validate that the annotated VCF is not empty
    if [ ! -s ${sample_id}.vep.annotated.vcf ]; then
        echo "Error: VEP annotation output is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}


process BCFTOOLS_MERGE {
	container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
	publishDir "${params.outdir}/BCFTOOLS_MERGE", mode: "copy"
	
    input:
    path vcfs

    output:
    path "merged_output.vcf.gz"
    path "merged_output.vcf.gz.tbi"
	path "variants_per_sample.tsv" 


    script:
    """
    vcfs_absolute=\$(for vcf in ${vcfs.join(' ')}; do realpath "\$vcf"; done)

    echo "Using absolute paths:" > debug_absolute_paths.log
    echo \$vcfs_absolute >> debug_absolute_paths.log

    bcftools merge \\
        \$vcfs_absolute \\
        -O z -o merged_output.vcf.gz
    tabix -p vcf merged_output.vcf.gz
	
	# Query merged VCF to generate a tabular summary
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' merged_output.vcf.gz > variants_per_sample.tsv
    """
}




process ANNOTATE_VARIANTS {
    tag "Annotate variants"
    
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf	// Filtered VCF file
	path index
	path tsv
    path snpEffJar  // Path to the SnpEff JAR file
    path snpEffConfig  // Path to the SnpEff configuration file
    path snpEffDbDir	// Path to the SnpEff database directory
	val genomedb

    output:
    path "annotated.vcf"
	path "annotated.summary.html"

    script:
	
    
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated.vcf

    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated.summary.html \
        ${vcf} > /dev/null
    """
}


process VCF_TO_TABLE {
    tag "Convert VCF to Table"

    container 'https://depot.galaxyproject.org/singularity/cyvcf2%3A0.31.1--py312h68a07e8_1'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf_file
	path html
	path script_file

    output:
    path "filtered_variants.csv"
    

    script:
    """
	
	
	
	python ${script_file} ${vcf_file}
	
    """
}

process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
	publishDir "${params.outdir}/annotations", mode: 'copy'
    input:
    path "input.vcf.gz"          // Input VCF file
    path "input.vcf.gz.tbi"      // Tabix index file for the VCF
	path tsv
    path "vep_cache"             // VEP cache directory
    path "clinvar.vcf.gz"        // ClinVar VCF file
    path "clinvar.vcf.gz.tbi"    // Tabix index file for ClinVar

    output:
    path "annotated_variants.vcf"          // Annotated VCF output
    path "annotated_variants.html"         // HTML summary report

    script:
    """
    vep --input_file input.vcf.gz \
        --output_file annotated_variants.vcf \
        --stats_file annotated_variants.html \
        --cache \
        --dir_cache vep_cache \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN
    """
}


process MULTIQC_REPORT {
    tag "MultiQC"
	container "https://depot.galaxyproject.org/singularity/multiqc%3A1.26--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_report_files 

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${qc_report_files.join(' ')} -o .
    """
}



workflow {
    // 1. Input channel for recalibrated BAM files
    aligned_bams = Channel.fromPath("/home/kothai/cq-git-sample/test/data/test/STAR/*.bam")
                           .map { bam -> tuple(bam.baseName.split('_Aligned')[0], bam) }
						   


	
	// Sort and index BAM files
    sorted_bams = SAMTOOLS_SORT_INDEX(aligned_bams)
	
	//Filter BAMS and index BAM files
	filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)
	
	// Generate alignment statistics
    alignment_stats = SAMTOOLS_FLAGSTAT(filtered_bams)
	
	
	
	// Mark duplicates
    marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
	
	// Split N CIGAR reads
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict) })
	
	//SAMTOOLS CALMD process
	calmd_ch = SAMTOOLS_CALMD(split_bams , params.genome, params.fasta_index )
	
	// Recalibrate and Apply BQSR in one step
    recalibrated_bams = GATK_RECALIBRATION(
        calmd_ch.map { tuple(it[0], it[1], it[2], params.genome, params.fasta_index, params.genome_dict, params.filtered_vcf, params.filtered_vcf_index) })
		
	// Convert BED to interval list
    interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist, params.genome, params.genome_dict)
	
	// Scatter the Interval List
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.genome_dict)
	
		
	//GATK HaplotypeCaller
	gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams, params.genome, params.fasta_index, params.genome_dict, scattered_intervals_ch)
				
	
	gvcf_output.view { "Raw GVCF output: $it" }
	
	//Step 15: provide stats
	bcftools_stats_ch = BCFTOOLS_STATS(gvcf_output)
	
	// Step 1: Filter individual VCF files
	filtered_individual_vcfs = GATK_VARIANT_FILTER(gvcf_output, params.genome, params.fasta_index, params.genome_dict)
	
	//Provide Stats
	filtered_vcf_stats = BCFTOOLS_QUERY(filtered_individual_vcfs) 

	// Step 2: Create a mapped collection of filtered VCF paths for merging
	filtered_vcf = filtered_individual_vcfs
                  .map { it[1] } // Extract paths to filtered VCF files
                  .collect()     // Collect them into a list

	filtered_vcf.view { "Filtered VCF output: $it" }

	// Step 3: Dynamically set paths based on the mode (test or actual)
	def snpEffJar = file(params.mode == 'test' ? './data/test/snpEff/snpEff.jar' : './data/actual/snpEff/snpEff.jar')
	def snpEffConfig = file(params.mode == 'test' ? './data/test/snpEff/snpEff.config' : './data/actual/snpEff.config')
	def snpEffDbDir = file(params.mode == 'test' ? './data/test/snpEff/snpEff/data' : './data/actual/snpEff/data')

	// Step 4: Conditional processing for merging or annotating individual files
	if (params.merge_vcf) {
		// Merge filtered VCFs
		merged_filtered_vcfs = BCFTOOLS_MERGE(filtered_vcf)

    // Annotate the merged VCF file snpeff
    annotated_merged_vcf = ANNOTATE_VARIANTS(merged_filtered_vcfs, snpEffJar, snpEffConfig, snpEffDbDir, params.genomedb)
	
	//Annotate the merged vcfs ensembl_vep
	annotated_merged_vcf_ensemblvep = ANNOTATEVARIANTS_VEP(merged_filtered_vcfs, params.vep_cache, params.clinvar, params.clinvartbi)

    println "Merging and annotating VCF files completed."
	
	// Create a table from the annotated merged VCF
    def vcf_to_table_script = file('convert_vcf_to_table.py') 
    table_creation = VCF_TO_TABLE(annotated_merged_vcf, vcf_to_table_script)

    println "Table creation from merged VCF completed."
	
	} else {
    
	// Annotate individual VCF files
    annotated_individual_vcfs = ANNOTATE_INDIVIDUAL_VARIANTS(filtered_individual_vcfs, snpEffJar, snpEffConfig, snpEffDbDir, params.genomedb)
	
	//Annotate individual VCF files ensemblvep
	annotated_individual_vcf_ensemblvep = ANNOTATE_INDIVIDUAL_VARIANTS_VEP (filtered_individual_vcfs, params.vep_cache, params.clinvar, params.clinvartbi)

    println "Individual VCF annotation completed."
    println "Table creation step skipped because merging is disabled."
}
	qc_outputs_ch = alignment_stats
    .map { it[1] } // Extract path from tuple
    .mix(
        bcftools_stats_ch.map { it[1] }, // Extract path
        filtered_vcf_stats.map { it[1] } // Extract path
    )
    .collect()


	multiqc_results = MULTIQC_REPORT(qc_outputs_ch)

	
	
	
	

}