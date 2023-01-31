process wcx_convert {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")

	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz")
	script:
		"""
		WisecondorX convert ${sampleID}/${sampleID}.bam ${sampleID}/${sampleID}.bam.wcx.npz
		"""
}

process wcx_predict {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz")
	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}_WCXpredict_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict_bins.bed"),file("${sampleID}/${sampleID}_WCXpredict.plots/*")

	script:
		"""
		WisecondorX --loglevel info predict ${sampleID}/${sampleID}.bam.wcx.npz ${params.wcx_reference} ${sampleID}/${sampleID}_WCXpredict --blacklist ${params.blacklist} --bed --plot --zscore ${params.wcxZscore} 
		"""
}

process wcx_predict_preface {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        //errorStrategy 'ignore'

        input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz")

        output:
		tuple val(sampleID),file("${sampleID}/${sampleID}_WCXpredict.preface_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict.preface_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_bins.bed")

        script:
                """
		WisecondorX --loglevel info predict ${sampleID}/${sampleID}.bam.wcx.npz ${params.preface_wcx_reference} ${sampleID}/${sampleID}_WCXpredict.preface --blacklist ${params.blacklist} --bed --zscore ${params.wcxZscore}
                """
}

process wcx_gender {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz")

	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}.wcx.gender.txt")

	script:
		"""
		WisecondorX gender ${sampleID}/${sampleID}.bam.wcx.npz ${params.wcx_reference} > ${sampleID}/${sampleID}.wcx.gender.txt
		"""
}

