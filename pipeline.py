import time
import os,sys
import subprocess
def oss(cmd):
	print cmd
	os.system(cmd)

root=sys.path[0]
genome="hg19"

if len(sys.argv)>2:
	genome=sys.argv[2]


genome_motifset=dict()
genome_motifset["hg19"]="tfh"
genome_motifset["mm9"]="tf"
genome_motifset["mm10"]="tfh"

motifdir="/broad/compbio/pouyak/motifs/verts/insts/"+genome+"/"+genome_motifset[genome]+"/Intergenic/optmm"
inputPeakfile=sys.argv[1]

genome_maskname=dict()
genome_maskname["hg19"]="268435455"
genome_maskname['mm9']="1048575"
genome_maskname["mm10"]="1099511627775"
maskname=genome_maskname[genome]


motifMatchFile=motifdir+"/8mer/"+maskname+"/matches/0.0.gz"
tmpdir="/broad/hptmp/zhizhuo/motifpipeline/"+os.path.basename(inputPeakfile).split(".")[0]
cmd="mkdir -p "+tmpdir
oss(cmd)
tmpfile=tmpdir+"/"+os.path.basename(inputPeakfile)

tmpfile_origwidth=tmpdir+"/"+os.path.basename(inputPeakfile)+".origwidth"
topN=10000
topNfilterStr="|awk '$6<"+str(topN)+"{print}'"
##reformat peak file
cmd="sh "+root+"/reformat.sh "+inputPeakfile+topNfilterStr+" > "+tmpfile
oss(cmd)
cmd="sh "+root+"/reformat_originalwidth.sh "+inputPeakfile+topNfilterStr+"  > "+tmpfile_origwidth
oss(cmd)


def IsGoodSize(fn):
	if not os.path.isfile(fn):
		return False
	statinfo = os.stat(fn)
	if statinfo.st_size>80000:
		return True
	else:
		return False 
	
##compute global enrichment score
gbdir=tmpfile+".global_enrichment_score"
GlobalEnrichmentFn=gbdir+"/enrichments.txt.gz"
if not IsGoodSize(GlobalEnrichmentFn):
	cmd="zsh "+root+"/Enricher.zzz.sh -d "+gbdir+" -o "+genome+" -i "+tmpfile_origwidth+" -C 0 -p  "+motifdir+" -P 0 -Z 1.96 -N 0 -n 0.0 -k "+maskname+"  -K 8"  ##set Number of Node =0  , make it local execution 
	oss(cmd+" &")

##grep motif-overlap
tmpfile2=tmpfile+".ol"
if not os.path.isfile(tmpfile2): 
	cmd="grep-overlap "+motifMatchFile+" "+tmpfile+"|grep -v 8mer_C > "+tmpfile2
	oss(cmd)

##compute positional and peak-rank score
PosRankEnrichmentFn=tmpfile2+".pickle.score"
if not os.path.isfile(PosRankEnrichmentFn): 
	cmd="python "+root+"/computePositionBias_RankBias.py "+tmpfile2
	oss(cmd)


##compute CombineAvgScore
round=0
while True:
	round+=1
	cmd="ps aux|grep zzhang|grep Enricher.zzz|wc -l"
	#if IsGoodSize(GlobalEnrichmentFn) or round>2000:
	output=subprocess.check_output(cmd, shell=True).strip()
	print output
	if output=="2":
		break
	time.sleep(3)
cmd="python "+root+"/combineAvgScore.py "+GlobalEnrichmentFn+" "+PosRankEnrichmentFn
oss(cmd)
#
inputSummaryPickleFn=PosRankEnrichmentFn+".summary.pickle"
cmd="python "+root+"/AcceptReject.py "+inputSummaryPickleFn
oss(cmd)

tfname_HGNCname=dict()
for line in open(root+"/data/name-mapping.txt"):
	comps=line.strip().split()
	tfname_HGNCname[comps[0].upper()]=comps[1]
TFname="unknownTF"
comps=os.path.basename(inputPeakfile).split("_")
if len(comps)>1:
	TFname=comps[1].upper()
	if TFname in tfname_HGNCname:
		TFname=tfname_HGNCname[TFname]
	TFname=TFname.replace("ALPHA","A").replace("BETA","B").replace("GAMMA","C").upper()
print TFname 
cmd="python "+root+"/makehtml.py "+tmpdir+" "+TFname+" > "+tmpdir+"/index.html"
oss(cmd)
