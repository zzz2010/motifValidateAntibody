import os,sys

def oss(cmd):
	print cmd
	os.system(cmd)

root=sys.path[0]
genome="hg19"
motifdir="/broad/compbio/pouyak/motifs/verts/insts/"+genome+"/tfh/Intergenic/optmm"
inputPeakfile="data/K562_ATF4_11815_Rep1_2_RepAll_peaks-2.bed"
maskname="268435455"
motifMatchFile=motifdir+"/8mer/"+maskname+"/matches/0.0.gz"
tmpdir="/tmp/motifpipeline/"+os.path.basename(inputPeakfile).split(".")[0]
cmd="mkdir -p "+tmpdir
oss(cmd)
tmpfile=tmpdir+"/"+os.path.basename(inputPeakfile)

##reformat peak file
cmd="sh "+root+"/reformat.sh "+inputPeakfile+"  > "+tmpfile
oss(cmd)

##compute global enrichment score
gbdir=tmpfile+".global_enrichment_score"
cmd="Enricher.sh -d "+gbdir+" -o "+genome+" -i "+tmpfile+" -C 0 -p  "+motifdir+" -P 0 -Z 1.96 -N 0 -n 0.0 -k "+maskname+"  -K 8"  ##set Number of Node =0  , make it local execution 
oss(cmd+" &")

##grep motif-overlap
tmpfile2=tmpfile+".ol"
if not os.path.isfile(tmpfile2): 
	cmd="grep-overlap "+motifMatchFile+" "+tmpfile+"|grep -v 8mer_C > "+tmpfile2
	oss(cmd)

##compute positional and peak-rank score
cmd="python "+root+"/computePositionBias_RankBias.py "+tmpfile2
oss(cmd)


##compute CombineAvgScore
GlobalEnrichmentFn=gbdir+"/enrichments.txt.gz"
PosRankEnrichmentFn=tmpfile2+".pickle.score"
cmd="python "+root+"/combineAvgScore.py "+GlobalEnrichmentFn+" "+PosRankEnrichmentFn
oss(cmd)

inputSummaryPickleFn=PosRankEnrichmentFn+".summary.pickle"
cmd="python "+root+"/AcceptReject.py "+inputSummaryPickleFn
oss(cmd)
