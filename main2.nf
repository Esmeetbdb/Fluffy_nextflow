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
}

if (params.skipline){
        Channel
                .fromPath( params.samplesheet )
                .splitCsv(header: true, skip: 1)
                .map{ row->row.Sample_ID }
                .unique()
                .set{ sample_ch }
}

//Load modules
include {  bwa_aln_R1; bwa_aln_R2; bwa_sampe; samtools_sort;picard_md} from './modules/preproccess.nf'
include {  fastqc_R1; fastqc_R2; collect_gc_bias; collect_insert_size; estimate_complexity} from './modules/qc.nf'

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
	collect_gcbias_bam_ch=collect_gc_bias(picard_md_output_ch)
	insert_size_output_ch=collect_insert_size(picard_md_output_ch)
	complexity_metrics_ch=estimate_complexity(picard_md_output_ch)


}

process wcx_convert {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") from wcx_convert_bam_ch

	output:
		set val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz") into wcx_convert_output_ch
	script:
		"""
		WisecondorX convert ${sampleID}/${sampleID}.bam ${sampleID}/${sampleID}.bam.wcx.npz
		"""
}

wcx_convert_output_ch.into { wcx_predict_input_ch; wcx_gender_input_ch; wcx_preface_predict_input_ch }

process wcx_predict {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz") from wcx_predict_input_ch
	output:
		set val(sampleID),file("${sampleID}/${sampleID}_WCXpredict_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict_bins.bed") into wcx_predict_output_ch

	script:
		"""
		WisecondorX --loglevel info predict ${sampleID}/${sampleID}.bam.wcx.npz ${params.wcx_reference} ${sampleID}/${sampleID}_WCXpredict --bed --zscore ${params.wcxZscore} 
		"""
}

process wcx_predict_preface {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz") from wcx_preface_predict_input_ch

        output:
		set val(sampleID),file("${sampleID}/${sampleID}_WCXpredict.preface_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict.preface_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_bins.bed") into wcx_preface_predict_output_ch

        script:
                """
		WisecondorX --loglevel info predict ${sampleID}/${sampleID}.bam.wcx.npz ${params.preface_wcx_reference} ${sampleID}/${sampleID}_WCXpredict.preface --bed --zscore ${params.wcxZscore}
                """
}

process preface_predict {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		set val(sampleID),file("${sampleID}/${sampleID}_WCXpredict.preface_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict.preface_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict.preface_bins.bed") from wcx_preface_predict_output_ch
	output:
		set val(sampleID),file("${sampleID}/${sampleID}_bins.bed.PREFACE.txt") into preface_predict_output_ch

	script:
		"""
		Rscript /bin/PREFACE-0.1.1/PREFACE.R predict --infile ${sampleID}/${sampleID}_WCXpredict.preface_bins.bed --model ${params.model} > ${sampleID}/${sampleID}_bins.bed.PREFACE.txt
		"""
}	
	
	
process wcx_gender {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam.wcx.npz") from wcx_gender_input_ch

	output:
		set val(sampleID),file("${sampleID}/${sampleID}.wcx.gender.txt") into wcx_gender_output_ch

	script:
		"""
		WisecondorX gender ${sampleID}/${sampleID}.bam.wcx.npz ${params.wcx_reference} > ${sampleID}/${sampleID}.wcx.gender.txt
		"""
}

process tiddit {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") from tiddit_bam_ch

        output:
                set val(sampleID),file("${sampleID}/${sampleID}.tiddit.tab") into tiddit_output_ch

        script:
                """
		tiddit --cov --bam ${sampleID}/${sampleID}.bam -z ${params.tiddit_binsize} -o ${sampleID}/${sampleID}.tiddit
		"""
}

process get_gctab {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        output:
                path("gc_tab") into gc_tab_output_ch

        script:
                """
		python ${params.script_folder}/AMYCNE/Generate_GC_tab.py --fa ${params.reference} --size ${params.tiddit_binsize} --n_mask > gc_tab
		"""
}

process run_amyce {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

        input:
		set val(sampleID),file("${sampleID}/${sampleID}.tiddit.tab") from tiddit_output_ch
                path("gc_tab") from gc_tab_output_ch


        output:
                set val(sampleID),file("${sampleID}/${sampleID}.tiddit.AMYCNE.tab") into AMYCNE_output_ch

        script:
                """
		python ${params.script_folder}/AMYCNE/AMYCNE.py --ff --coverage ${sampleID}/${sampleID}.tiddit.tab --gc gc_tab --Q ${params.Q} --scaling ${params.scaling} --intercept ${params.intercept} >${sampleID}/${sampleID}.tiddit.AMYCNE.tab
		"""
}

summarize_input_ch = Channel.fromPath( "${params.output}" )
multiQC_input_ch= Channel.fromPath( "${params.output}" )



summarize_temp_1_ch = AMYCNE_output_ch.cross(wcx_gender_output_ch).map{
	it ->  [it[0][0],it[0][1],it[1][1]]
    }
summarize_temp_2_ch = summarize_temp_1_ch.cross(wcx_predict_output_ch).map{
	it -> [it[0][0],it[0][1],it[0][2],it[1][1],it[1][2],it[1][3],it[1][4]]
    }


process summarize {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path (out_dir) from summarize_input_ch
		set val(sampleID), file("${sampleID}/${sampleID}.tiddit.AMYCNE.tab"), file("${sampleID}/${sampleID}.wcx.gender.txt"),file("${sampleID}/${sampleID}_WCXpredict_aberrations.bed"),file("${sampleID}/${sampleID}_WCXpredict_chr_statistics.txt"),file("${sampleID}/${sampleID}_WCXpredict_segments.bed"),file("${sampleID}/${sampleID}_WCXpredict_bins.bed") from summarize_temp_2_ch
	output:
		file("test_output.csv") into final_output_ch

	script:
		"""
		python ${params.script_folder}/generate_csv.py --folder ${out_dir} --samplesheet ${params.samplesheet} --Zscore ${params.Zscore} --minCNV ${params.minCNV} --maxGCD ${params.maxGCD} --maxATD ${params.maxATD} --maxbin2bin ${params.maxbin2bin} --maxdup ${params.maxdup} --minreads ${params.minreads} > test_output.csv
		"""
}

multiQC_temp_1_ch = insert_size_output_ch.cross(complexity_metrics_ch).map{
	it -> [it[0][0],it[0][1],it[0][2],it[1][1]]
    }
multiQC_temp_2_ch = multiQC_temp_1_ch.cross(gc_bias_output_ch).map{
	it -> [it[0][0],it[0][1],it[0][2],it[0][3],it[1][1],it[1][2],it[1][3]]
    }

process multiQC {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'

	input:
		path(out_dir) from multiQC_input_ch
		set val(sampleID),file("${sampleID}/${sampleID}.insert_size_metrics.txt"),file("${sampleID}/${sampleID}.insert_size_histogram.pdf"),file("${sampleID}/${sampleID}.complex_metrics.txt"),file("${sampleID}/${sampleID}.gc_bias_metrics.txt"),file("${sampleID}/${sampleID}.gc_bias_metrics.pdf"),file("${sampleID}/${sampleID}.gc.summary.tab") from multiQC_temp_2_ch

	output:
		file("${out_dir}/multiqc_report.html") into multiQC_output_ch

	script:
		"""
		multiqc ${out_dir}/**/ --outdir $out_dir
		"""
}
