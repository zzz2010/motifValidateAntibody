peakfile=$1

sort -k5gr $peakfile|awk   'BEGIN{OFS="\t";}{if(NF<3||$3<$2){$3=$2+1}$2=int(($2+$3)/2);print ".",$1,$2-100,$2+100,$5,NR}'|rtool sort  

