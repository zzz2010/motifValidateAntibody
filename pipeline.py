import time
import os,sys

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

motifdir="/broad/compbio/pouyak/motifs/verts/insts/"+genome+"/"+genome_motifset[genome]+"/Intergenic/optmm"
inputPeakfile=sys.argv[1]

genome_maskname=dict()
genome_maskname["hg19"]="268435455"
genome_maskname['mm9']="1048575"
maskname=genome_maskname[genome]


motifMatchFile=motifdir+"/8mer/"+maskname+"/matches/0.0.gz"
tmpdir="/broad/hptmp/zhizhuo/motifpipeline/"+os.path.basename(inputPeakfile).split(".")[0]
cmd="mkdir -p "+tmpdir
oss(cmd)
tmpfile=tmpdir+"/"+os.path.basename(inputPeakfile)

##reformat peak file
cmd="sh "+root+"/reformat.sh "+inputPeakfile+"  > "+tmpfile
oss(cmd)

def IsGoodSize(fn):
	if not os.path.isfile(fn):
		return False
	statinfo = os.stat(fn)
	print statinfo.st_size
	if statinfo.st_size>30000:
		return True
	else:
		return False 
	
##compute global enrichment score
gbdir=tmpfile+".global_enrichment_score"
GlobalEnrichmentFn=gbdir+"/enrichments.txt.gz"
if not IsGoodSize(GlobalEnrichmentFn):
	cmd="zsh "+root+"/Enricher.zzz.sh -d "+gbdir+" -o "+genome+" -i "+tmpfile+" -C 0 -p  "+motifdir+" -P 0 -Z 1.96 -N 0 -n 0.0 -k "+maskname+"  -K 8"  ##set Number of Node =0  , make it local execution 
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
	if IsGoodSize(GlobalEnrichmentFn) or round>2000:
		break
	time.sleep(3)
cmd="python "+root+"/combineAvgScore.py "+GlobalEnrichmentFn+" "+PosRankEnrichmentFn
oss(cmd)
#
inputSummaryPickleFn=PosRankEnrichmentFn+".summary.pickle"
cmd="python "+root+"/AcceptReject.py "+inputSummaryPickleFn
oss(cmd)

TFname="unknownTF"
comps=os.path.basename(inputPeakfile).split("_")
if len(comps)>1:
	TFname=comps[1]
print TFname 
cmd="python "+root+"/makehtml.py "+tmpdir+" "+TFname+" > "+tmpdir+"/index.html"
oss(cmd)
