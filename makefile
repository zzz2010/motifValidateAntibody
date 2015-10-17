data/name-mapping.txt:~/compbio/projects/centdist2/result/2015-02-12/CentdistScoreForPouyaPipelineOutput/name-mapping.txt
	cut -f 2,3 $<|sed 's/_/\t/1'|cut -f 1,3|sed 's/_/\t/1'|cut -f 1,2|sort -u > $@
