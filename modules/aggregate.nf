process summarize {	
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path (out_dir)
		val (fluffy_results)
		file (samplesheet)
	output:
		file("*NIPT.csv")

	script:
		"""
		python ${params.script_folder}/generate_csv.py --folder ${out_dir} --samplesheet ${samplesheet} --Zscore ${params.Zscore} --minCNV ${params.minCNV} --maxGCD ${params.maxGCD} --maxATD ${params.maxATD} --maxbin2bin ${params.maxbin2bin} --maxdup ${params.maxdup} --minreads ${params.minreads}
		"""
}

process multiQC {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path(out_dir)
		val(qc_files)

	output:
		file("${out_dir}/multiqc_report.html")

	script:
		"""
		multiqc ${out_dir}/**/ --outdir $out_dir
		"""
}
