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
		file("multiqc_report.html")

	script:
		"""
		multiqc ${out_dir}/**/
		"""
}

process make_deliverables_yaml {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path(out_dir)
		val(multiqc)
		val(summary)

	output:
		file("deliverables.yaml")

	script:
		"""
		python ${params.script_folder}/make_deliverables_yaml.py ${out_dir}
		"""
}
