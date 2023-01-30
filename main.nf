nextflow.enable.dsl=2

params.samplesheet = false
params.fastq = false
params.output = false
params.skipline = false
params.prefix = "Sample_"
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

//Load modules
include {bwa_aln_R1; bwa_aln_R2; bwa_sampe; samtools_sort;picard_md} from './modules/preproccess.nf'
include {fastqc_R1; fastqc_R2; collect_gc_bias; collect_insert_size; estimate_complexity} from './modules/qc.nf'
include {wcx_convert; wcx_predict; wcx_predict_preface; wcx_gender} from './modules/wisecondor.nf'
include {preface_predict; tiddit; get_gctab; run_amyce} from './modules/fetalfraction.nf'
include {summarize; multiQC} from './modules/aggregate.nf'

//main workflow
workflow{

	//Align and sort
	bwa_aln_R1_ch=bwa_aln_R1(sample_ch)
	bwa_aln_R2_ch=bwa_aln_R2(sample_ch)

	fastqc_R1_output_ch=fastqc_R1(sample_ch)
	fastqc_R2_output_ch=fastqc_R2(sample_ch)

	bwa_sampe_input_ch = bwa_aln_R1_ch.cross(bwa_aln_R2_ch).map{
        	it ->  [it[0][0],it[0][1],it[1][1]]
	}

	bwa_sampe_output_ch=bwa_sampe(bwa_sampe_input_ch)

	samtools_sort_output_ch=samtools_sort(bwa_sampe_output_ch)

	//picard QC
	picard_md_output_ch=picard_md(samtools_sort_output_ch)
	gc_bias_output_ch=collect_gc_bias(picard_md_output_ch)
	insert_size_output_ch=collect_insert_size(picard_md_output_ch)
	complexity_metrics_ch=estimate_complexity(picard_md_output_ch)

	//run wisecondorX modules
	wcx_convert_output_ch=wcx_convert(picard_md_output_ch)
	wcx_predict_output_ch=wcx_predict(wcx_convert_output_ch)
	wcx_gender_output_ch=wcx_gender(wcx_convert_output_ch)
	wcx_preface_predict_output_ch=wcx_predict_preface(wcx_convert_output_ch)


	//run fetal fraction prediction tools
	gc_tab_output_ch=get_gctab()
	preface_predict_output_ch=preface_predict(wcx_preface_predict_output_ch)
	tiddit_output_ch=tiddit(picard_md_output_ch)	
	AMYCNE_output_ch=run_amyce(tiddit_output_ch,gc_tab_output_ch)

	//aggregate results
	summarize_input_ch = Channel.fromPath( "${params.output}" )
	multiQC_input_ch= Channel.fromPath( "${params.output}" )

	summarize_temp_1_ch = AMYCNE_output_ch.cross(wcx_gender_output_ch).map{
		it ->  [it[0][0],it[0][1],it[1][1]]
	}
	summarize_temp_2_ch = summarize_temp_1_ch.cross(wcx_predict_output_ch).map{
		it -> [it[0][0],it[0][1],it[0][2],it[1][1],it[1][2],it[1][3],it[1][4]]
	}

	multiQC_temp_1_ch = insert_size_output_ch.cross(complexity_metrics_ch).map{
		it -> [it[0][0],it[0][1],it[0][2],it[1][1]]
	}

	multiQC_temp_2_ch = multiQC_temp_1_ch.cross(gc_bias_output_ch).map{
		it -> [it[0][0],it[0][1],it[0][2],it[0][3],it[1][1],it[1][2],it[1][3]]
	}

	summary_output_ch=summarize(summarize_input_ch,summarize_temp_2_ch)
	multiQC_output_ch=multiQC(multiQC_input_ch,multiQC_temp_2_ch)

}
