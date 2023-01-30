process summarize {	
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path (out_dir)
		tuple val(sampleID), file("${sampleID}/${sampleID}.tiddit.AMYCNE.tab"), file("${sampleID}/${sampleID}.wcx.gender.txt"),file("${sampleID}/${sampleID}_WCXpredict_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict_bins.bed")
		file (samplesheet)
	output:
		file("test_output.csv")

	script:
		"""
		python ${params.script_folder}/generate_csv.py --folder ${out_dir} --samplesheet ${samplesheet} --Zscore ${params.Zscore} --minCNV ${params.minCNV} --maxGCD ${params.maxGCD} --maxATD ${params.maxATD} --maxbin2bin ${params.maxbin2bin} --maxdup ${params.maxdup} --minreads ${params.minreads} > test_output.csv
		"""
}

process multiQC {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path(out_dir)
		path(qc_files)

	output:
		file("${out_dir}/multiqc_report.html")

	script:
		"""
		multiqc ${out_dir}/**/ --outdir $out_dir
		"""
}
