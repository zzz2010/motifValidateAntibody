import os,sys
import glob
import numpy as np
from sklearn import metrics
import scipy.stats
import gzip

GlobalEnrichmentFn=sys.argv[1]
PosRankEnrichmentFn=sys.argv[2]

motifNames_clust=dict()
TFNameSet=set()
gene2Protein=dict()
for line in open("data/gene2Protein.txt"):
        comps=line.strip().split()
        gene2Protein[comps[0].upper()]=comps[1].upper()

		
motifsWinstance=set()
for line in open("data/uniqueMotif.txt"):
        motifsWinstance.add(line.strip().upper())

cid=-1
cid_name=dict()
multigroup_TF=set()
for line in open("data/motifs-clust-names.txt"):
        cid+=1
        comps=line.strip().replace(";","\t").split()
        motifset=set(comps)
        tfset=set()
        for m in motifset:
                m=m.upper()
                tfname=m.split("_")[0]
                if m not in motifsWinstance:
                        continue
                if tfname in motifNames_clust and motifNames_clust[tfname]!=cid: ##mean there is exist cluster contain it
                        #print tfname
                        multigroup_TF.add(tfname)
                if tfname not in motifNames_clust:
                        motifNames_clust[tfname]=set()
                motifNames_clust[tfname].add(cid)
                motifNames_clust[m]=cid
                if tfname in TFNameSet or len(tfset)<30:
                        tfset.add(tfname)
        if len(tfset)==0:
                continue
        clustName="|".join(list(tfset))
        cid_name[cid]=clustName+"|"+str(cid)

		
TotalMotifCandidates=set()

for line in open("data/uniqueMotif.txt"):
        m=line.strip().upper()
        tfname=m.split("_")[0]
        if len(tfname)<2:
                continue
        if m not in motifNames_clust:
                continue
        if len(cid_name[motifNames_clust[m]])<2:
                continue
        TotalMotifCandidates.add(motifNames_clust[m])

		
def getRankandLabel(fn):
	mncid_bestmotif=dict()  ##record tf_motifcid to individual motif mapping
        lines=open(fn).readlines()
        #print fn
        numScoreTypes=len(lines[0].split())-1
        nrow=max(len(lines),len(TotalMotifCandidates))
        score_mat=np.zeros((nrow,numScoreTypes))
        Y=[0]*nrow
        i=0
        motifNameList=list()
        visited=dict()
        for line in lines:
                comps=line.replace("nan","-10").strip().split()
                motifname=comps[0].upper()
                if motifname not in motifNames_clust:
                        continue
                motifname_cid=motifNames_clust[motifname]
                if len(cid_name[motifname_cid])<2: ##filter motif with name only one character
                        continue
		mncid=motifname.split("_")[0]+"_"+str(motifname_cid)
                if motifname_cid  in visited:
                        ii=visited[motifname_cid]
			oldscore=0
			newscore=0
                        for j in range(0,numScoreTypes):
				oldscore+=score_mat[ii,j]
				newscore+=float(comps[j+1])
                                score_mat[ii,j]=max(score_mat[ii,j] ,float(comps[j+1]))
			if mncid not in mncid_bestmotif:
				mncid_bestmotif[mncid]=motifname
			elif oldscore<newscore:
				mncid_bestmotif[mncid]=motifname
                        continue
                else:
			mncid_bestmotif[mncid]=motifname
                        motifNameList.append(motifname_cid)
                        visited[motifname_cid]=i
                        for j in range(0,numScoreTypes):
                                score_mat[i,j]=float(comps[j+1])
                i+=1
        aa=np.amin(score_mat[0:i,:],axis=0)
        for m in TotalMotifCandidates:
                if m in visited:
                        continue
                motifNameList.append(m)
                score_mat[i,:]=aa
                i+=1

        rank_mat=np.zeros((i,numScoreTypes))
        Zscore_mat=np.zeros((i,numScoreTypes))
        for j in range(0,numScoreTypes):
                rank_mat[:,j]=scipy.stats.rankdata(score_mat[0:i,j])
                Zscore_mat[:,j]=scipy.stats.mstats.zscore(score_mat[0:i,j])

  ##make sure only count the best match motif group in each dataset
        ranklist=rank_mat[0:i,0].copy()
        Y=np.array(Y[0:i])
        ranklist[np.argwhere(Y==0)]=0
        Y=np.zeros(i)
        Y[np.argmax(ranklist)]=1
        Y=list(Y)

        ##
        return (rank_mat,Y[0:i],motifNameList,Zscore_mat[0:i,:],mncid_bestmotif)

####load globalenrichment score########
dataset_motifscore=dict()
for line in gzip.open(GlobalEnrichmentFn):
        comps=line.strip().split()
        if comps[2]=="+":
                if comps[3]!=".":
                        continue
                d=comps[3]
                m=comps[1]
                if m not in motifNames_clust:
                        continue
                cid=motifNames_clust[m]
                score=float(comps[13])
                if d not in dataset_motifscore:
                        dataset_motifscore[d]=dict()
                if cid not in  dataset_motifscore[d]:
                        dataset_motifscore[d][cid]=float(score)
                else:
                        dataset_motifscore[d][cid]=max(float(score),dataset_motifscore[d][cid])

						
import cPickle as pickle
combined_rank_mat=None
combined_zscore_mat=None
combined_Y=None
combined_chipseqmotif=list()
outf=open(PosRankEnrichmentFn+".summary.txt","w")
chipseqmotif_tf=list()
chipseq_tf=list()
fl=PosRankEnrichmentFn
bname=os.path.basename(fl).replace(".pickle","").replace(".score","").replace(".bed.ol","")
TFname="."
X,Y,motifnamelist,score_mat,mncid_bestmotif=getRankandLabel(fl)
##expand 2 columns
X=np.column_stack((X,np.zeros((X.shape[0], 2)))) #ranking matrix
score_mat=np.column_stack((score_mat,np.zeros((score_mat.shape[0], 2)))) #score matrix
gbcol=2
avgcombin_col=gbcol+1
l=-1
fail=0
for m in motifnamelist:
	l+=1
	if m in dataset_motifscore['.']:
	        score_mat[l,gbcol]=dataset_motifscore['.'][m]
	else:   
	        fail+=1
chipseq_tf.append(TFname)
X[:,gbcol]=scipy.stats.rankdata(score_mat[:,gbcol])
score_mat[:,gbcol]=scipy.stats.mstats.zscore(score_mat[:,gbcol])
score_mat[:,avgcombin_col]=scipy.stats.mstats.zscore(score_mat[:,gbcol]+score_mat[:,0]+score_mat[:,1] )
X[:,avgcombin_col]=scipy.stats.rankdata(score_mat[:,avgcombin_col])
for i in range(len(motifnamelist)):
        oline=bname+"\t"+cid_name[motifnamelist[i]]+"\t"+"\t".join(map(str,score_mat[i,:]))
        outf.write(oline+"\n")
        chipseqmotif_tf.append(TFname)
        combined_chipseqmotif.append(bname+"\t"+cid_name[motifnamelist[i]])
if combined_rank_mat is None:
        combined_rank_mat=X
        combined_Y=Y
        combined_zscore_mat=score_mat
else:
        combined_rank_mat=np.vstack((combined_rank_mat,X))
        combined_zscore_mat=np.vstack((combined_zscore_mat,score_mat))
        combined_Y+=Y
outf.close()


pickle.dump( [combined_Y,combined_rank_mat,combined_chipseqmotif,combined_zscore_mat,mncid_bestmotif], open( PosRankEnrichmentFn+".summary.pickle", "wb" ) )	
