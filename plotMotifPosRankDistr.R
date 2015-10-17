options(stringsAsFactors=FALSE)

Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
inputOlFile=Args[2]
motif=Args[3]
outdir=Args[4]
cmd=sprintf("grep %s %s|tr '|' '\t'|tr ' '  '\t' ",motif,inputOlFile)
df=read.table(pipe(cmd),sep="\t")

system(sprintf("mkdir -p %s",outdir))

mpos=as.integer((df[,3]+df[,4])/2)
ppos=as.integer((df[,10]+df[,11])/2)

distances=abs(mpos-ppos)
ranks=df[,13]

pdf(sprintf("%s/%s.pdf",outdir,motif),height=4,width=8)
par(mfrow=c(1,2))
hist(distances,10,col='red',xlab="Distance to Peak Center(bp)" ,main="Position Distribution" )
hist(ranks,10,col='blue',xlab="Peak Rank",main="Peak-Rank Distribution" )

dev.off()















