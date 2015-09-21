import os,sys
import pickle
overlapfile=sys.argv[1]
peak_name=os.path.basename(overlapfile)
print >> sys.stderr,  peak_name
maxRange=1000
maxRange2=0
maxRank=0

PositionCount=dict()
RankCount=dict()
expName=os.path.splitext(os.path.dirname(overlapfile))[1]
fn=overlapfile+".pickle"
if True: #not os.path.isfile(fn):
	for line in open(overlapfile):
		comps=line.strip().split()
		motifName=comps[0].replace("_8mer","")
		motif_pos=(int(comps[2])+int(comps[3]))/2
		peakcenter=(int(comps[8])+int(comps[9]))/2
		distance=int(abs(motif_pos-peakcenter))
#		if distance> maxRange-1:
#			continue
		peakwidth=(int(comps[9])-int(comps[8]))/2
		if peakwidth>maxRange2:
			maxRange2=peakwidth
		peakRank=int(comps[11].split('|')[0])
		if maxRank< peakRank:
			maxRank=peakRank
		if motifName not in PositionCount:
			PositionCount[motifName]=dict()
			RankCount[motifName]=set()
		if peakRank not in RankCount[motifName]:
			RankCount[motifName].add(peakRank)
		if distance not in PositionCount[motifName]:
			PositionCount[motifName][distance]=0
		PositionCount[motifName][distance]+=1
	
	pickle.dump([PositionCount,RankCount,maxRank,maxRange2], open(fn, "wb" ) )
else:
	[PositionCount,RankCount,maxRank,maxRange2]=pickle.load(open(fn, "rb" ))

print >> sys.stderr,  maxRank,maxRange2
import scipy.stats
import math
def RankBiasScore_wilcox(rankCount1,maxRank1):
	countArr=list()
	uniformRank=list()
	for i in range(maxRank1):
		if i in rankCount1:
			countArr.append(i)
		else:
			uniformRank.append(i)
	
	res=scipy.stats.ranksums(countArr,uniformRank)
	a=-res[0]
	if math.isnan(a):
		a=-2
	return a



def PositionalBiasScore_binom(PositionCount1,maxRange1):
	inCount=0
	outCount=1
	winsize=50
	for i in PositionCount1:
		if i< winsize:
			inCount+=PositionCount1[i]
		else:
			outCount+=PositionCount1[i]
	p=float(winsize)/maxRange1
	n=(inCount+outCount)
	z=(inCount-n*p)/math.sqrt(n*p*(1-p))
	return z #-scipy.stats.binom.logsf(inCount,inCount+outCount,float(winsize)/maxRange1)


outf=open(fn+".score",'w')
for motif in RankCount:
	pScore=PositionalBiasScore_binom(PositionCount[motif],maxRange2)
	rScore2=RankBiasScore_wilcox(RankCount[motif],maxRank)
	outf.write(motif+"\t"+str(pScore)+"\t"+str(rScore2)+"\n")
outf.close()
