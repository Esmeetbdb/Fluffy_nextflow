process fastqc_R1 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	val(sampleID)
        output:
               	tuple val(sampleID),file("${sampleID}/*fastqc.*")
        script:
               	"""
                mkdir -p ${sampleID}
                zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R1*fastq.gz | fastqc stdin:${sampleID}_R1 -o ${sampleID}
                """
}

process fastqc_R2 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	val(sampleID)
        output:
               	tuple val(sampleID), file("${sampleID}/*fastqc.*")

        script:
               	"""
                mkdir -p ${sampleID}
                zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R2*fastq.gz | fastqc stdin:${sampleID}_R2 -o ${sampleID}
                """
}

process collect_gc_bias {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")

	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.gc_bias_metrics.txt"),file("${sampleID}/${sampleID}.gc_bias_metrics.pdf"),file("${sampleID}/${sampleID}.gc.summary.tab")

	script:
		"""
		picard CollectGcBiasMetrics I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.gc_bias_metrics.txt CHART=${sampleID}/${sampleID}.gc_bias_metrics.pdf S=${sampleID}/${sampleID}.gc.summary.tab R=${params.reference} VALIDATION_STRINGENCY=LENIENT TMP_DIR=${params.tmpdir}
		"""
}

process collect_insert_size {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")
	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.insert_size_metrics.txt"),file("${sampleID}/${sampleID}.insert_size_histogram.pdf")

	script:
		"""
		picard CollectInsertSizeMetrics I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.insert_size_metrics.txt H=${sampleID}/${sampleID}.insert_size_histogram.pdf VALIDATION_STRINGENCY=LENIENT M=0.5 TMP_DIR=${params.tmpdir}
		"""
}

process estimate_complexity {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")
	
	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.complex_metrics.txt")

	script:
		"""
		picard EstimateLibraryComplexity I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.complex_metrics.txt VALIDATION_STRINGENCY=LENIENT TMP_DIR=${params.tmpdir}
		"""
}
