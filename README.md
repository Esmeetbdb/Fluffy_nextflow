# Fluffy
Fluffy is an NIPT analysis pipeline written in Python and Nextflow. Fluffy predicts chromosomal trisomies using fetal fraction estimation by AMYCNE. Fluffy also computes multiQC data of the input sequence files. Fluffy outputs one csv file with a line per input sample which summarizes the main results per sample, including the Z-score of chromosomes prone to trisomies.

# Run
```
Fluffy can be directly run in Nextflow using:
    nextflow run Fluffy.nf --samplesheet /path/to/samplesheet.csv --fastq /path/to/fastq/folder --output /path/to/output/folder -c config.conf
```
If the Samplesheet file has a line before the header add the following to the command: 
```
--skipline true
```

### Input
The Fluffy pipeline requires two distinct inputs: the samplesheet and the a folder with fastq files per sample.
##### Samplesheet

The samplesheet should have a header with the following column names:
```
FCID,Lane,Sample_ID,SampleRef,index,index2,SampleName,Control,Recipe,Operator,Sample_Project
```
There should be one comm-separated line for each sample to be analyzed

##### Fastq folder
The folder with fastq files should contain one folder for each sample to be analysed. The folder names should be 'Sample_{Sample_ID}'.
Inside this folder should be all the fastq files for the sample. There should be at least one file for the forward reads and one file for reverse reads. If multiple lanes are used during sequencing this is not a problem. The file names should be '\*{sample_ID}\*\_R1\*fastq.gz' for forward reads and '\*{sample_ID}\*\_R2\*fastq.gz' for reverse reads.
**All fastq files should be gzipped**

### Output
Fluffy produces several output files. 
- summary.csv This file contains a header and one line per sample summarizing all the data extracted from the sample during analysis
- multiqc.html This file contains the multiQC report of all samples
- One folder per sample with folder name {Sample_ID}. In this folder there will be the output file of each separate step in the Fluffy pipeline

# Installation
Dependencies:
```
Nextflow
```
Install the dependencies then download Fluffy.
```
git clone --recursive https://github.com/Esmeetbdb/Fluffy_nextflow
```
Next, download the singularity images specified in the config.conf file.
```
singularity pull singularity pull docker://link/to/singularity/image
```
For the Fluffy singularity image use.
```
singularity pull library://jeisfeldt/default/fluffy:sha256.dbef92cd5eab8558c2729f73a191d73a7576a24e9bb44dde7372c0cd405c4ef6
```
The path to all singularity containers must be specified in the config file.

You will need to download/create the following files:
```
Reference fasta (indexed using bwa)
WisecondorX reference files (created using the reference mode)
PREFACE model file (optional)
```

# Config file
An example config file can be found in the configs folder. Default parameters are specied. Read the comments in the config.conf file for more information on the individual variables.
