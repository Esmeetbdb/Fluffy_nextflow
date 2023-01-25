process bwa_aln_R1 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	val(sampleID)

        output:
               	tuple val(sampleID),file("${sampleID}/${sampleID}_R1.sai")

        script:
               	def R1 = "<( zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R1*fastq.gz )"

                """
                mkdir -p ${sampleID}
                bwa aln -n 0 -k 0 -t ${task.cpus} ${params.reference} ${R1} > ${sampleID}/${sampleID}_R1.sai
                """
}

process bwa_aln_R2 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	val(sampleID)

        output:
                tuple val(sampleID),file("${sampleID}/${sampleID}_R2.sai")

        script:
               	def R2 = "<( zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R2*fastq.gz )"

                """
                mkdir -p ${sampleID}
                bwa aln -n 0 -k 0 -t ${task.cpus} ${params.reference} ${R2} > ${sampleID}/${sampleID}_R2.sai
                """
}

process bwa_sampe {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	tuple val(sampleID),file("${sampleID}/${sampleID}_R1.sai"),file("${sampleID}/${sampleID}_R2.sai")
        output:
               	tuple val(sampleID),file("${sampleID}/${sampleID}_sampe.sam")

        script:
               	def R1 = "<( zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R1*fastq.gz )"
                def R2 = "<( zcat ${params.fastq}/${params.prefix}${sampleID}/*${sampleID}*_R2*fastq.gz )"

                """
                bwa sampe -n -1 ${params.reference} ${sampleID}/${sampleID}_R1.sai ${sampleID}/${sampleID}_R2.sai ${R1} ${R2} > ${sampleID}/${sampleID}_sampe.sam
                """
}


process samtools_sort {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
              	tuple val(sampleID),file("${sampleID}/${sampleID}_sampe.sam")

        output:
               	tuple val(sampleID),file("${sampleID}/${sampleID}.tmp.bam")

        script:
               	"""
                samtools sort --output-fmt bam ${sampleID}/${sampleID}_sampe.sam -@ 8 -m 2G -T ${sampleID}/${sampleID}.tmp > ${sampleID}/${sampleID}.tmp.bam
                """
}

process picard_md {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID) ,file("${sampleID}/${sampleID}.tmp.bam")

	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")

	script:
		"""
		picard MarkDuplicates I=${sampleID}/${sampleID}.tmp.bam O=${sampleID}/${sampleID}.bam M=${sampleID}/${sampleID}.md.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT  TMP_DIR=${params.tmpdir}
		rm ${sampleID}/${sampleID}.tmp.bam
		"""
}
