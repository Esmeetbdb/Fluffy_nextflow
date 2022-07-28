params.samplesheet = false
params.fastq = false
params.output = false

if (!params.samplesheet || !params.fastq || !params.output){
	println "missing required parameters"
	println "main.nf --samplesheet <samplesheet_csv> --fastq <fastq_folder> --output <output_folder>"
	exit o
}

Channel
	.fromPath( params.samplesheet )
	.splitCsv(header: true, skip: 1)
	.map{ row-> row.Sample_ID }
	.unique()
	.set{ sample_ch }

sample_ch.into { sample_ch_R1; sample_ch_r2; fastqc_R2_input_ch; fastqc_R1_input_ch }

process bwa_aln_R1 {
	publishDir "${params.output}", mode: 'copy', overwrite: true
	errorStrategy 'ignore'
        cpus 1

	input:
		val(sampleID) from sample_ch_R1

	output:
		set val(sampleID),file("${sampleID}/${sampleID}_R1.sai") into bwa_aln_R1_ch

	script:
		def R1 = "<( zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R1*fastq.gz )"

		"""
		mkdir -p ${sampleID}
		bwa aln -n 0 -k 0 -t ${task.cpus} ${params.reference} ${R1} > ${sampleID}/${sampleID}_R1.sai
		"""
}

process fastqc_R1 {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

        input:
		val(sampleID) from fastqc_R1_input_ch
	output:
		set val(sampleID),file("${sampleID}/*fastqc.*") into fastqc_R1_output_ch
	script:
		"""
		mkdir -p ${sampleID}
		zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R1*fastq.gz | fastqc stdin:${sampleID}_R1 -o ${sampleID}
		"""
}		

process bwa_aln_R2 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
	errorStrategy 'ignore'
        cpus 1

        input:
                val(sampleID) from sample_ch_r2

        output:
                set val(sampleID),file("${sampleID}/${sampleID}_R2.sai") into bwa_aln_R2_ch

        script:
                def R2 = "<( zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R2*fastq.gz )"

                """
		mkdir -p ${sampleID}
                bwa aln -n 0 -k 0 -t ${task.cpus} ${params.reference} ${R2} > ${sampleID}/${sampleID}_R2.sai
                """
}

process fastqc_R2 {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

        input:
                val(sampleID) from fastqc_R2_input_ch
	output:
		set val(sampleID), file("${sampleID}/*fastqc.*") into fastqc_R2_output_ch

        script:
                """
		mkdir -p ${sampleID}
		zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R2*fastq.gz | fastqc stdin:${sampleID}_R2 -o ${sampleID}
                """
}

bwa_sampe_input_ch = bwa_aln_R1_ch.cross(bwa_aln_R2_ch).map{
        it ->  [it[0][0],it[0][1],it[1][1]]
    }


process bwa_sampe {
        publishDir "${params.output}", mode: 'copy', overwrite: true
	errorStrategy 'ignore'
        cpus 8

        input:
		set val(sampleID),file("${sampleID}/${sampleID}_R1.sai"),file("${sampleID}/${sampleID}_R2.sai") from bwa_sampe_input_ch
        output:
                set val(sampleID),file("${sampleID}/${sampleID}_sampe.sam") into bwa_sampe_output_ch

        script:
                def R1 = "<( zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R1*fastq.gz )"
                def R2 = "<( zcat ${params.fastq}/Sample_${sampleID}/*${sampleID}*_R2*fastq.gz )"

                """
                bwa sampe -n -1 ${params.reference} ${sampleID}/${sampleID}_R1.sai ${sampleID}/${sampleID}_R2.sai ${R1} ${R2} > ${sampleID}/${sampleID}_sampe.sam
		"""
}


process samtools_sort {
        publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 8

        input:
                set val(sampleID),file("${sampleID}/${sampleID}_sampe.sam") from bwa_sampe_output_ch

        output:
                set val(sampleID),file("${sampleID}/${sampleID}.tmp.bam") into samtools_sort_output_ch

        script:
                """
                samtools sort --output-fmt bam ${sampleID}/${sampleID}_sampe.sam -@ 8 -m 2G -T ${sampleID}/${sampleID}.tmp > ${sampleID}/${sampleID}.tmp.bam
                """
}

process picard_md {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

	input:
		set val(sampleID) ,file("${sampleID}/${sampleID}.tmp.bam") from samtools_sort_output_ch

	output:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") into picard_md_output_ch

	script:
		"""
		picard MarkDuplicates I=${sampleID}/${sampleID}.tmp.bam O=${sampleID}/${sampleID}.bam M=${sampleID}/${sampleID}.md.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
		rm ${sampleID}/${sampleID}.tmp.bam
		"""
}

picard_md_output_ch.into { wcx_convert_bam_ch; tiddit_bam_ch; collect_gcbias_bam_ch;  collect_insertsize_bam_ch; collect_estimatecomplexity_bam_ch }

process collect_gc_bias {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

	input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") from collect_gcbias_bam_ch

	output:
		set val(sampleID),file("${sampleID}/${sampleID}.gc_bias_metrics.txt"),file("${sampleID}/${sampleID}.gc_bias_metrics.pdf"),file("${sampleID}/${sampleID}.gc.summary.tab") into gc_bias_output_ch

	script:
		"""
		picard CollectGcBiasMetrics I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.gc_bias_metrics.txt CHART=${sampleID}/${sampleID}.gc_bias_metrics.pdf S=${sampleID}/${sampleID}.gc.summary.tab R=${params.reference} VALIDATION_STRINGENCY=LENIENT TMP_DIR=${params.tmpdir}
		"""
}

process collect_insert_size {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

        input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") from collect_insertsize_bam_ch
	output:
		set val(sampleID),file("${sampleID}/${sampleID}.insert_size_metrics.txt"),file("${sampleID}/${sampleID}.insert_size_histogram.pdf") into insert_size_output_ch

	script:
		"""
		picard CollectInsertSizeMetrics I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.insert_size_metrics.txt H=${sampleID}/${sampleID}.insert_size_histogram.pdf VALIDATION_STRINGENCY=LENIENT M=0.5 TMP_DIR=${params.tmpdir}
		"""
}

process estimate_complexity {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

        input:
		set val(sampleID),file("${sampleID}/${sampleID}.bam"),file("${sampleID}/${sampleID}.bai"),file("${sampleID}/${sampleID}.md.txt") from collect_estimatecomplexity_bam_ch
	
	output:
		set val(sampleID),file("${sampleID}/${sampleID}.complex_metrics.txt") into complexity_metrics_ch

	script:
		"""
		picard EstimateLibraryComplexity I=${sampleID}/${sampleID}.bam O=${sampleID}/${sampleID}.complex_metrics.txt VALIDATION_STRINGENCY=LENIENT TMP_DIR=${params.tmpdir}
		"""
}

process wcx_convert {
	publishDir "${params.output}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        cpus 1

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
        cpus 1

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
        cpus 1

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
        cpus 1

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
        scratch true
        cpus 1

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
        cpus 1

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
        cpus 1

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
        cpus 1

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
        cpus 1

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
        cpus 1

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
