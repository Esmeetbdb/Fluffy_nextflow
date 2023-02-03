import sys

coverages={}
total_bins=0.0
first=True
chromosomes=set(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"])
coverage_list=list( map(lambda x: x/100.0, range(0, 100, 1)) )
m=1
coverages[m]=0
for c in coverage_list:
	coverages[c]=0

for line in open(sys.argv[1]):
	content=line.strip().split()

	if first:
		first=False
		continue

	if float(content[3]) == 0 or float(content[4]) < 10:
		continue
	if not content[0].replace("chr","") in chromosomes:
		continue

	total_bins+=1
	v=round(float(content[3]),2)
	if v >=m:
		coverages[m]+=1 

	else:
		coverages[v]+=1



print("# id: Tiddit coverage summary")
print("# section_name: \'Tiddit coverage summary\'")
print("# description: \'This output is described in the file header. Any MultiQC installation will understand it without prior configuration.\'")
print("# format: \'tsv\'")
print("# plot_type: \'linegraph\'")
print("# pconfig:")
print("#    id: 'custom_bargraph_w_header'")
print("#    title: Per bin average coverage distribution")
print("#    ylab: 'Fraction of bins'")
print("#    height: 400")


for c in sorted(coverages.keys()):
	print("{}\t{}".format(c,coverages[c]/total_bins))
