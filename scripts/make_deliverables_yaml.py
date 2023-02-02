import sys
import yaml
import glob

#python make_deliverables.yaml <output_dir>

output_dir=sys.argv[1]
deliverables={}
deliverables["files"]=[]

summary=glob.glob(f"{output_dir}/*.csv", recursive=False)[0]
project_name=summary.split(".csv")[0].split("/")[-1]

multiqc_path=glob.glob(f"{output_dir}/*.html", recursive=False)[0]

deliverables["files"].append({"format":"csv", "id":project_name,"path":str(summary),"step":"summarise_batch","tag":"NIPT_csv"})
deliverables["files"].append({"format":"html", "id":project_name,"path":str(multiqc_path),"step":"summarise_batch","tag":"MultiQC"})

samples=glob.glob(f"{output_dir}/**/", recursive=False)
for sample in samples:

	sample_id=sample.strip("/").split("/")[-1]

	abberations=glob.glob(f"{sample}/*_WCXpredict_aberrations.bed", recursive=False)
	if abberations:
		deliverables["files"].append({"format":"bed", "id":sample_id,"path":str(abberations[0]),"step":"wisecondorX","tag":"Wisecondor_aberrations"})


	chr_stat=glob.glob(f"{sample}/*WCXpredict_statistics.txt", recursive=False)
	if chr_stat:
		deliverables["files"].append({"format":"txt", "id":sample_id,"path":str(chr_stat[0]),"step":"wisecondorX","tag":"Wisecondor_chromosome_statistics"})

	plots=glob.glob(f"{sample}/*_WCXpredict.plots", recursive=False)
	if plots:
		deliverables["files"].append({"format":"folder", "id":sample_id,"path":str(plots[0]),"step":"wisecondorX","tag":"Wisecondor_plots"})

	tiddit_cov=glob.glob(f"{sample}/*.tiddit.tab", recursive=False)
	if tiddit_cov:
		deliverables["files"].append({"format":"bed", "id":sample_id,"path":str(tiddit_cov[0]),"step":"tiddit","tag":"tiddit_coverage"})

f=open("deliverables.yaml","w")
f.write(yaml.dump(deliverables))
f.close()
