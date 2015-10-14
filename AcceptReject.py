import math
import gzip
import os,sys
import cPickle as pickle
import numpy as np
from scipy.stats import nbinom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=sys.path[0]
maxRank=283
inputSummaryPickleFn=sys.argv[1]
modelParamsTuple=pickle.load(open(root+"/cistrome_model.pickle",'rb'))
modelparameter_z1=modelParamsTuple[0]
modelparameter_z0=modelParamsTuple[1]
motif_support=modelParamsTuple[2]

##draw decision figure##
def NBPara2Arr(param,n):
        #print param
        Arr=np.zeros((n,1))
        for i in range(int(n)):
                Arr[i]=nbinom.pmf( i,1, 1-param) #change failure prob
        Arr*=1/np.sum(Arr)
        return Arr
def getMotifLogOdd(z1Param,z0Param,maxRank):
        z1Arr=NBPara2Arr(z1Param,maxRank)
        z0Arr=NBPara2Arr(z0Param,maxRank)
        logORarr=np.zeros((len(z0Arr),1))
        for i in range(len(z0Arr)):
                logORarr[i]=np.log(z1Arr[i])-np.log(z0Arr[i])
        return logORarr

def getKLdivergence(z1Param,z0Param,maxRank):
        z1Arr=NBPara2Arr(z1Param,maxRank)
        z0Arr=NBPara2Arr(z0Param,maxRank)
        return scipy.stats.entropy(z1Arr,z0Arr)

def getMaxlogOdd(vec):
        return np.amax(np.absolute(vec))

def simulationOddsRatio(z1Param,z0Param,maxRank,n_simulation,isPositive):
        z1Arr=np.ravel(NBPara2Arr(z1Param,maxRank))
        z0Arr=np.ravel(NBPara2Arr(z0Param,maxRank))
        if isPositive:
                ranks=np.random.choice(maxRank,n_simulation,replace=True,p=z1Arr )
        else:   
                ranks=np.random.choice(maxRank,n_simulation,replace=True,p=z0Arr )
        return np.array([np.log(z1Arr[r]/z0Arr[r]) for r in ranks])

scoreCol=3
def odd2prob(x):
        r=math.exp(x)
        return r/(r+1)
result=pickle.load(open( inputSummaryPickleFn, "rb" ) )

def Arr2str(a):
	return "\t".join(map(str,a))

badMotifs=set()
for line in open(root+"/data/badRankMotif.txt"):
	badMotifs.add(line.strip())	

i=0
Y=result[0]
ranks=result[1]
chipseq_motif=result[2]
zscores=result[3]
mncid_bestmotif=result[4]

dataset_logOdd=dict()
cm_strs=chipseq_motif[0].split()
chip=cm_strs[0].split("_")[0]
chipdataset=cm_strs[0]
motifs=cm_strs[1].split("|")
motif_cid=motifs[len(motifs)-1]
outf=open(inputSummaryPickleFn+".motifresult",'w')
header="Motif_cid\tIsBadMotif\tPositionBiasZscore\tPeakRankBiasZscore\tGlobalEnrichmentZscore\tCombineAvgScore\tEnrichmentRank(total:"+str(maxRank)+")\tbestMotif"
outf.write(header+"\n")
tfLogOdd=dict()
tfBestScore=dict()
tfBestMotif=dict()
tfBestMotifScore=dict()
for j in range(len(chipseq_motif)-1):
	score_vec=zscores[j]
	rank_vec=ranks[j]
        mc_str=chipseq_motif[j].split("\t")[1]
	toks=mc_str.split("|") ##motiflist + cid
	cid=toks[len(toks)-1]
	for motif in toks[0:(len(toks)-1)]:
		m=motif+"_"+cid
		isbad="Y" if m in badMotifs else "N"
		tf=motif
		if m not in mncid_bestmotif:
			continue
		bestmotif=mncid_bestmotif[m]
	        r=maxRank-rank_vec[scoreCol]
		outf.write(m+"\t"+isbad+"\t"+Arr2str(score_vec)+"\t"+str(r+1)+"\t"+bestmotif+"\n")
	        if m not in modelparameter_z0:
	                continue
	        if m in modelparameter_z1:
	                z1Param=modelparameter_z1[m]
	                support=motif_support[m]
	        else:
	                z1Param=modelparameter_z1["other"]
		if motif=="CTCF":
			print r, z1Param,modelparameter_z0[m]
	        logOdds=getMotifLogOdd(z1Param,modelparameter_z0[m],maxRank)[r]
		if tf not in tfLogOdd:
			tfLogOdd[tf]=float(logOdds)
			tfBestScore[tf]=float(logOdds)
			tfBestMotif[tf]=m
			tfBestMotifScore[tf]=Arr2str(score_vec)
		else:
			tfLogOdd[tf]+=float(logOdds)
			if tfBestScore[tf]<float(logOdds):
				tfBestScore[tf]=float(logOdds)
				tfBestMotif[tf]=m
				tfBestMotifScore[tf]=Arr2str(score_vec)
outf.close()

tf_validatability=dict()
i=-1
for line in open(root+"/data/TF_validatability.txt"):
	i+=1
	if i==0:
		continue
	comps=line.strip().split()
	tf_validatability[comps[0]]=comps[1]+"\t"+comps[2]

outf=open(inputSummaryPickleFn+".tfresult",'w')
header="TFName\tAcceptLogOdds\tAcceptProb\tBestMotif\tPositionBiasZscore\tPeakRankBiasZscore\tGlobalEnrichmentZscore\tCombineAvgScore\tRejectCapability\tAcceptCapability"
outf.write(header+"\n")
for tf in tfBestScore:
	logOdds=tfLogOdd[tf]
	if tf not in tf_validatability: #that mean not in training data
		continue
	outf.write(Arr2str([tf,float(logOdds),odd2prob(logOdds)])+"\t"+tfBestMotif[tf]+"\t"+(tfBestMotifScore[tf])+"\t"+tf_validatability[tf]+"\n")
outf.close()

