nextflow.enable.dsl=2

params.samplesheet = false
params.fastq = false
params.output = false
params.skipline = false
params.prefix = ""
params.version = "1.0.0"

println "running Fluffy version ${params.version}"

if (!params.samplesheet || !params.fastq || !params.output){
	println "missing required parameters"
	println "main.nf --samplesheet <samplesheet_csv> --fastq <fastq_folder> --output <output_folder>"
	exit 1
}

//Read samplesheet
if (!params.skipline){
	Channel
		.fromPath( params.samplesheet )
		.splitCsv(header: true)
		.map{ row->row.Sample_ID }
		.unique()
		.set{ sample_ch }
}else{
	Channel
		.fromPath( params.samplesheet )
		.splitCsv(header: true, skip: 1)
		.map{ row->row.Sample_ID }
		.unique()
		.set{ sample_ch }
}


if(params.output[0] ==  "/"){
	output_absolute_path=params.output
}else{
	output_absolute_path="${launchDir}/${params.output}"
}
println output_absolute_path

//Load modules
include {bwa_aln_R1; bwa_aln_R2; bwa_sampe; samtools_sort;picard_md} from './modules/preproccess.nf'
include {fastqc_R1; fastqc_R2; collect_gc_bias; collect_insert_size; estimate_complexity; samtools_flagstat; samtools_stat; coverage_summary} from './modules/qc.nf'
include {wcx_convert; wcx_predict; wcx_predict_preface; wcx_gender} from './modules/wisecondor.nf'
include {preface_predict; tiddit; get_gctab; run_amyce} from './modules/fetalfraction.nf'
include {summarize; multiQC; make_deliverables_yaml} from './modules/aggregate.nf'

samplesheet=file(params.samplesheet)

//main workflow
workflow{
	//Preprocess:align,sort and mark duplicates
	bwa_aln_R1_ch = bwa_aln_R1(sample_ch)
	bwa_aln_R2_ch = bwa_aln_R2(sample_ch)

	bwa_sampe_input_ch = bwa_aln_R1_ch.join(bwa_aln_R2_ch)
	bwa_sampe_output_ch = bwa_sampe(bwa_sampe_input_ch)

	samtools_sort_output_ch = samtools_sort(bwa_sampe_output_ch)
	picard_md_output_ch = picard_md(samtools_sort_output_ch)


	//run wisecondorX modules
	wcx_convert_output_ch = wcx_convert(picard_md_output_ch)
	wcx_predict_output_ch = wcx_predict(wcx_convert_output_ch)
	wcx_gender_output_ch = wcx_gender(wcx_convert_output_ch)
	wcx_preface_predict_output_ch = wcx_predict_preface(wcx_convert_output_ch)

	//run fetal fraction prediction tools
	gc_tab_output_ch = get_gctab()
	preface_predict_output_ch = preface_predict(wcx_preface_predict_output_ch)
	tiddit_output_ch = tiddit(picard_md_output_ch)	
	AMYCNE_output_ch = run_amyce(tiddit_output_ch,gc_tab_output_ch)

	//QC
	fastqc_R1_output_ch = fastqc_R1(sample_ch)
	fastqc_R2_output_ch = fastqc_R2(sample_ch)

	gc_bias_output_ch = collect_gc_bias(picard_md_output_ch)
	insert_size_output_ch = collect_insert_size(picard_md_output_ch)
	complexity_metrics_ch = estimate_complexity(picard_md_output_ch)

	samtools_flagstat_ch = samtools_flagstat(picard_md_output_ch )
	samtools_stat_ch = samtools_stat(picard_md_output_ch )

	coverage_summary_ch=coverage_summary(tiddit_output_ch)

	//aggregate results
	output_path = Channel.fromPath( "${params.output}" )

	summarize_ch = AMYCNE_output_ch.join(wcx_gender_output_ch)
	summarize_ch = summarize_ch.join(wcx_predict_output_ch)
	summarize_ch = summarize_ch.join(preface_predict_output_ch)
	summarize_ch = summarize_ch.join(insert_size_output_ch)
	summarize_ch = summarize_ch.join(gc_bias_output_ch)
	summarize_ch = summarize_ch.collect{it[1..-1]}

	multiQC_ch = picard_md_output_ch.join(gc_bias_output_ch)
	multiQC_ch = multiQC_ch.join(insert_size_output_ch)
	multiQC_ch = multiQC_ch.join(complexity_metrics_ch)
	multiQC_ch = multiQC_ch.join(fastqc_R1_output_ch)
	multiQC_ch = multiQC_ch.join(fastqc_R2_output_ch)
	multiQC_ch = multiQC_ch.join(samtools_flagstat_ch)
	multiQC_ch = multiQC_ch.join(samtools_stat_ch)
	multiQC_ch = multiQC_ch.join(coverage_summary_ch)
	multiQC_ch = multiQC_ch.collect{it[1..-1]}

	summary_output_ch = summarize(output_path,summarize_ch,samplesheet)
	multiQC_output_ch = multiQC(output_path,multiQC_ch)

	deliverables_yaml_ch = make_deliverables_yaml( output_absolute_path ,multiQC_output_ch,summary_output_ch)

}
