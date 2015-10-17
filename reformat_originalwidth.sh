peakfile=$1

sort -k5gr $peakfile|awk   'BEGIN{OFS="\t";}{if(NF<3||$3<$2){$3=$2+1};print ".",$1,$2,$3,$5,NR}'|rtool sort  

