process preface_predict {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		tuple val(sampleID),file("${sampleID}/${sampleID}_WCXpredict.preface_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict.preface_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_bins.bed")
	output:
		tuple val(sampleID),file("${sampleID}/${sampleID}_bins.bed.PREFACE.txt")

	script:
		"""
		Rscript ${params.script_folder}/PREFACE/PREFACE.R predict --infile ${sampleID}/${sampleID}_WCXpredict.preface_bins.bed --model ${params.model} > ${sampleID}/${sampleID}_bins.bed.PREFACE.txt
		"""
}	
	
process tiddit {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt")

        output:
                tuple val(sampleID),file("${sampleID}/${sampleID}.tiddit.tab")

        script:
                """
		tiddit --cov --bam ${sampleID}/${sampleID}.bam -z ${params.tiddit_binsize} -o ${sampleID}/${sampleID}.tiddit
		"""
}

process get_gctab {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        output:
                path("gc_tab")

        script:
                """
		python ${params.script_folder}/AMYCNE/Generate_GC_tab.py --fa ${params.reference} --size ${params.tiddit_binsize} --n_mask > gc_tab
		"""
}

process run_amyce {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		tuple val(sampleID),file("${sampleID}/${sampleID}.tiddit.tab")
                path("gc_tab")


        output:
                tuple val(sampleID),file("${sampleID}/${sampleID}.tiddit.AMYCNE.tab")

        script:
                """
		python ${params.script_folder}/AMYCNE/AMYCNE.py --ff --coverage ${sampleID}/${sampleID}.tiddit.tab --gc gc_tab --Q ${params.Q} --scaling ${params.scaling} --intercept ${params.intercept} >${sampleID}/${sampleID}.tiddit.AMYCNE.tab
		"""
}
