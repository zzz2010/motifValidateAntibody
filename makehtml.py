from decimal import *
import os,sys
import glob
import operator


root=sys.path[0]
dataPath=root+"/data/"
dataUrl="http://www.broadinstitute.org/~zzhang/motifpipeline/data/"

def oss(cmd):
	os.system(cmd+" >/dev/null 2>&1")
resultDir=sys.argv[1] #"/broad/hptmp/zhizhuo/motifpipeline/K562_ATF4_11815_Rep1_2_RepAll_peaks-2/" #sys.argv[1]
selTF=sys.argv[2]  #"ATF4"
tfresultFn=glob.glob(resultDir+"/*tfresult")[0]
motifresultFn=glob.glob(resultDir+"/*motifresult")[0]

motifclust_bestmotif=dict()
motifclust_rank=dict()

i=-1
for line in open(motifresultFn):
	i+=1
	if i==0:
		continue
	comps=line.strip().split()
	motifclust_bestmotif[comps[0]]=comps[len(comps)-1]	
	motifclust_rank[comps[0]]=comps[len(comps)-2]

def writeTag(tag,content):
        return "<"+tag+">"+content+"</"+tag+">\n"

def imageHtml(image_src,title):
        ret_str=image_src.replace("/home/unix/zzhang/public_html","http://www.broadinstitute.org/~zzhang").replace("/broad/hptmp/zhizhuo/motifpipeline","http://www.broadinstitute.org/~zzhang/motifpipeline/jobdata")
        ret_str="<div  style=\"width:218px; overflow-x:hidden;\"><embed src="+ret_str+" width=218px height=200px  onload=lzld(this)></embed>"+title+"</div>\n"
        return ret_str	

def list2TDs(str_list):
	ret=""
	for s in str_list:
		ret+=writeTag("td",s)
	return ret

def list2THs(str_list):
        ret=""
        for s in str_list:
                ret+=writeTag("th",s)
        return ret


def motifPosRankPlot(bestmotif):
	olFile=glob.glob(resultDir+"/*.ol")[0]
	outfn=resultDir+"/plots/"+bestmotif+".pdf"
	if not os.path.isfile(outfn):
		cmd="Rscript "+root+"/plotMotifPosRankDistr.R "+olFile+" "+bestmotif+" "+resultDir+"/plots/"	
		oss(cmd)
	ret_str=outfn.replace("/home/unix/zzhang/public_html","http://www.broadinstitute.org/~zzhang").replace("/broad/hptmp/zhizhuo/motifpipeline","http://www.broadinstitute.org/~zzhang/motifpipeline/jobdata")
	ret_str="<div  style=\"width:420px; overflow-x:hidden;\"><embed src="+ret_str+" width=420px  height=200px  onload=lzld(this)></embed> Motif Position and  Peak-Rank Distribution </div>\n"
	return ret_str 

def MotifClustRow(motifclust):
	outstr=""
        bestmotif=motifclust_bestmotif[motifclust]
	clustname=writeTag("td","cluster:<br>"+motifclust)
	imageFn=dataPath+"/motifDist/"+motifclust+".pdf"
	if os.path.isfile(imageFn):
		negbinomModel="<td colspan=2>"+imageHtml(imageFn.replace(dataPath,dataUrl),"Trained Model for Enrichment Rank")+"</td>"
	else:
		negbinomModel="<td colspan=2>"+"No Training Data"+"</td>"
        enrichRank=writeTag("td","Enrichment Rank:<br>"+motifclust_rank[motifclust])
        motiflogo="<td colspan=1>"+motifLogo(bestmotif)+"</td>"
        posrank_plot="<td colspan=4>"+ motifPosRankPlot(bestmotif)+"</td>"
        outstr=writeTag("tr",clustname+enrichRank+negbinomModel+motiflogo+posrank_plot )
	return outstr

def TFImageTRs(tf):
	imageFiles=glob.glob(dataPath+"/motifDist/"+tf+"*.pdf")
	outstr=""
	for imageFn in imageFiles:
		motifclust=os.path.basename(imageFn).replace(".pdf","")		
		outstr+=MotifClustRow(motifclust)
	return outstr

def motifLogo(motif):
	logourl="http://www.broadinstitute.org/~zzhang/logos_pdf/"+motif+".big.pdf"
        ret_str="<div  style=\"width:218px; overflow-x:hidden;\"><embed  src="+logourl+" width=218px  height=200px  onload=lzld(this)></embed>best motif: "+motif+"</div>\n"
	return ret_str

def TFTable(lines):
	header="TFName  AcceptLogOdds   AcceptProb      BestMotifClusterId       PositionBiasZscore      PeakRankBiasZscore      GlobalEnrichmentZscore  CombineAvgScore RejectCapability        AcceptCapability"
	selFields=[0,1,2,8,9,3,7,4,5,6]
	headerlist=header.split()
	outstr=list2THs([headerlist[s] for s in selFields])
	outstr=writeTag("tr",outstr) ##first line
	for line in lines:
		comps=line.strip().split()
		tdstr=list2TDs([comps[s] for s in selFields])
		outstr+=writeTag("tr",tdstr) 
		outstr+=TFImageTRs(comps[0])		
	return "<table border=1>"+outstr+"</table>"

def TFTable_noimage(lines):
        header="TFName  AcceptLogOdds   AcceptProb      BestMotifClusterId       PositionBiasZscore      PeakRankBiasZscore      GlobalEnrichmentZscore  CombineAvgScore RejectCapability        AcceptCapability"
        selFields=[0,1,2,8,9,3,7,4,5,6]
        headerlist=header.split()
        outstr=list2THs([headerlist[s] for s in selFields])
        outstr=writeTag("tr",outstr) ##first line
        for line in lines:
                comps=line.strip().split()
                tdstr=list2TDs([comps[s] for s in selFields])
                outstr+=writeTag("tr",tdstr)
        return "<table border=1>"+outstr+"</table>"

def MotifClustTable(lines):
	header="MotifClusterId       IsWrongMotif      PositionBiasZscore      PeakRankBiasZscore      GlobalEnrichmentZscore  CombineAvgScore EnrichmentRank(total:283)       bestMotif"
	selFields=[0,1,2,3,4,5]
	headerlist=header.split()
	outstr=list2THs([headerlist[s] for s in selFields])
	for line in lines:
		comps=line.strip().split()
		tdstr=list2TDs([comps[s] for s in selFields])
		outstr+=writeTag("tr",tdstr)
		outstr+=MotifClustRow(comps[0])
	return "<table border=1>"+outstr+"</table>"
	

htmlContent=""
CSSstr=" <style type=\"text/css\">img[alt=small] { float: center; width: 300px;}img[alt=large]{float: center; width: 800px;}</style>"
htmlContent+=CSSstr

if resultDir.endswith("/"):
        resultDir=resultDir[:-1]
Title=os.path.basename(resultDir)
htmlContent+=writeTag("head",writeTag("title",Title))

htmlContent+="<body>"
htmlContent+=" <script src=\"http://www.broadinstitute.org/~zzhang/motifpipeline/data/lazyload.min.js\"></script>"
#####make summary page #####
## selected TF statistic and bestMotif ## 
htmlContent+="<hr>"
htmlContent+=writeTag("H2","Selected TF: "+selTF)
selTFlines=list()
tf_accprb=dict()
acceptTFs=set()
rejectTFs=set()
decision_cutoff=0.8
i=-1
for line in open(tfresultFn):
	i+=1
	if i==0:
		continue
	comps=line.strip().split()
	if selTF in comps[0]:
		selTFlines.append(line)
	tf_accprb[comps[0]]=float(comps[2])
	prob=float(comps[2])
	if prob>decision_cutoff:
		acceptTFs.add(comps[0])
	elif prob<(1-decision_cutoff):
		rejectTFs.add(comps[0])
if len(selTFlines)>0:
	htmlContent+=TFTable(selTFlines)
if len(selTFlines)==0: ##possible bad motif
	i=-1
	sellines=list()
	for line in open(motifresultFn):
		i+=1
		if i==0:
			continue
		comps=line.strip().split()
		if selTF == comps[0].split("_")[0]:
			sellines.append(line)
	if len(sellines)>0:
		htmlContent+=MotifClustTable(sellines)
##accept TF list
htmlContent+="<hr>"
htmlContent+=writeTag("H2","Accepted TFs")

selTFlines=list()
for line in open(tfresultFn):
	comps=line.strip().split()
	if comps[0] in acceptTFs:
        	selTFlines.append(line)
htmlContent+=TFTable(selTFlines)

##top 10 best MotifClust ##
htmlContent+="<hr>"
htmlContent+=writeTag("H2","Top Enriched Motifs")
selMotifClust=set()
###top motif clust not metioned before ###
for mc in motifclust_rank:
	tf=mc.split("_")[0]
	if tf in acceptTFs or selTF in tf:
		continue
	if float(motifclust_rank[mc])<=10:
		selMotifClust.add(mc)

selMotiflines=dict()
for line in open(motifresultFn):
	comps=line.strip().split()
	if comps[0] in selMotifClust:
		selMotiflines[line]=float(comps[len(comps)-2])
import operator
sorted_lines = sorted(selMotiflines.items(), key=operator.itemgetter(1))
htmlContent+=MotifClustTable([x[0] for x in sorted_lines])

##reject TF list
htmlContent+="<hr>"
htmlContent+=writeTag("H2","Rejected TFs")
selTFlines=list()
for line in open(tfresultFn):
        comps=line.strip().split()
        if comps[0] in rejectTFs:
                selTFlines.append(line)
htmlContent+=TFTable_noimage(selTFlines)


htmlContent+="</body>"
print htmlContent

