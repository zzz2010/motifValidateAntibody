#!/bin/zsh

# intended to eventually do most analysis that would be interesting all automatically
# and can easily be browsed

# add new types of plots as they are used

HSTR=0
BACKREGIONS=Intergenic
CROPBREGIONS=0
ZCORR=1.5
FDRFILTER=0.05
ONPLOTS=3
MATCHESNAMES=$(seq -w 0.0 0.1 1.01 | tr '\n' ' ')
ENRICHBACKS=""
KEEPWITHMAX=0
LOWMODE=0
PRODUCERCF=1
NUMNODES=500
CONFIDZ=1

#TODO: add support for non-pwm motifs?
#TODO: add options to merge regions/motifs

if [ $# -eq 0 ]; then
		echo >&2 "USAGE: $0 [OPTIONS] < DIR/regions.txt
 -d     Output directory [required]
 -o     Organism. See GetOrgVars.sh for details. [required]
 -r     List of regions to consider for the motif figures (filters stdin) [default: (all)]
 -S     Do all analysis only on the first strand (0/1) [default: $HSTR]
 -i     Input file of regions (blank for stdin) [default: $INREGIONS]
 -C     Produce region based counts [default: $PRODUCERCF]
 -N     Number of nodes to use [default: $NUMNODES]

 -p     optmm directory to use [either -p or -m required]
 -m     Motif file to use [if used, then -n, -k -K are required; matching done with first of -R]
 -l     Low disk mode; skips production of matches (forces -n 0.0 -k 0 -C 0) (0/1) [default: $LOWMODE]

 -n     Names of the files to work on (must be gzipped) (blank to infer) [default: $MATCHESNAMES]
 -k     Match masks to use [default: infer automatically]
 -K     Kmer cutoffs to use (Z for supplied values in -m) [default: infer automatically]
 -X     Match motifs at most stringent cutoff when too short for -K [default: $KEEPWITHMAX]
 -R     Background regions (Name|[file]|[name_filter] where default file is ORGREGIONDIR/clean/Name.txt). Names should not be repeated, cannot include | [default: $BACKREGIONS]
 -Z     Z value used to correct confidence measurement [default: $CONFIDZ]

 -b     List of foreground regions to use as background for enrichment (in addition to + and ++ for everything and union of foreground); NOTE: must already be superset of regions considered! 

 -z     z-cutoff used for computing enrichments [default: $ZCORR]

 -P     Plots to produce (bitstring) [default: $ONPLOTS]

for -P 1=2^0: one plot per cutoff, showing all motifs/experiments
 -a     FDR filter -a value [default: $FDRFILTER]

for -P 2=2^1: one plot per motif, showing all cutoffs/experiments
 -M     List of motifs to consider [default: (all)]

for -P 4=2^2: cumulative plot indicating regions recovered by a matching motif (requires -C 1)
for -P 8=2^3: one plot showing all matching motifs for each factor
 -A     Alias file whose first semi-colon separated first column should be treated as aliases"
	exit 1
fi

while getopts d:R:o:p:n:k:K:S:M:z:r:a:m:P:i:A:b:X:l:C:N:Z: o
do      case "$o" in
		d)		DIR="$OPTARG";;
		i)		INREGIONS="$OPTARG";;
		R)		BACKREGIONS="$OPTARG";;
		o)		ORG="$OPTARG";;
		p)		OPTMMDIR="$OPTARG";;
		n)		MATCHESNAMES="$OPTARG";;
		k)		MATCHMASKS="$OPTARG";;
		m)		MOTIFFILE="$OPTARG";;
		M)		MOTIFLIST="$OPTARG";;
		r)		REGIONLIST="$OPTARG";;
		K)		KMERS="$OPTARG";;
		S)		HSTR="$OPTARG";;
		z)		ZCORR="$OPTARG";;
		a)		FDRFILTER="$OPTARG";;
		P)		ONPLOTS="$OPTARG";;
		A)		ALIASFILE="$OPTARG";;
		b)		ENRICHBACKS="$OPTARG";;
		X)		KEEPWITHMAX="$OPTARG";;
		l)		LOWMODE="$OPTARG";;
		C)		PRODUCERCF="$OPTARG";;
		N)		NUMNODES="$OPTARG";;
		Z)		CONFIDZ="$OPTARG";;
        esac
done

if [ $LOWMODE = 1 ]; then
	MATCHESNAMES=0.0
	PRODUCERCF=0
	MATCHMASKS=0
fi

if [ -z "$DIR" -o -z "$OPTMMDIR" -o "$LOWMODE" = 1 ] && [ -z "$DIR" -o -z "$ORG" -o -z "$MOTIFFILE" -o -z "$MATCHMASKS" -o -z "$KMERS" ]; then
	echo >&2 "ERROR: either -d -p -l 0 or -d -o -m -n -k -K are required."; exit 1
fi

if [ ! -z "$OPTMMDIR" -a ! -z "$MOTIFFILE" ]; then	
	echo >&2 "ERROR: cannot have both -p and -m."; exit 1
fi

if [ "$LOWMODE" = 1 -a $BACKREGIONS != ${BACKREGIONS/ */} ]; then
	# this restriction could be lifted but in order to guarantee the same output 
	# all the background regions would have to be filtered against the first one
	echo >&2 "ERROR: only one -R permitted with -l 1."; exit 1
fi 

mkdir -p $DIR

########################### PROCESS INPUT REGIONS
# process regions file creating a special one that is the union of all the regions
# add unique identifier to the end of each line so we couple duplicates accurately
if [ ! -f $DIR/regions.txt.gz ]; then 
       echo "making $DIR/regions.txt.gz" 
	if [ ! -z "$INREGIONS" ]; then
		gunzip -cf $INREGIONS
	else
		cat
	fi | awk -vOFS="\t" -vRL="$REGIONLIST" '
		BEGIN{
			split(RL, A, / /)
			
			for (i in A) R[A[i]]
		}
		
		(RL != "" && !($1 in R)){next}
		
		{
			print $1, $2, $3, $4, $5 != "" ? $5 : ".", ++X
			print "++", $2, $3, $4, $5 != "" ? $5 : ".", ++X
		}
		' | gzip -c > $DIR/regions.txt.gz
	 echo "done: $DIR/regions.txt.gz"
fi

if [ -z "$REGIONLIST" ]; then
	# use sort to make REGIONLIST ordered
	REGIONLIST=$(gunzip -cf $DIR/regions.txt | awk '!($1 in X){print $1; X[$1]}' | sort | awk '$0!="++"{L=L " " $0}; END{print substr(L,2)}')
fi


########################### SETUP REGIONS FILES FOR FUTURE COMPARE ANALYSIS
function rname () {
	echo ${1/|*/}
}

for r in $(echo $BACKREGIONS); do 
	mkdir -p $DIR/compare/$(rname $r)
	if [ ! -f $DIR/compare/$(rname $r)/bregions.txt.gz ]; then 
		echo "making $DIR/compare/$(rname $r)/bregions.txt.gz"
		# $r has format Name|[file]|[name_filter]
		if [ -z "$(echo $r | pcut -t '|' -f 2)" ]; then
			t=$(GetOrgVars.sh $ORG ORGREGIONDIR)/clean/$(rname $r).txt
		else
			t=$(echo $r | pcut -t '|' -f 2)		
		fi

		if [ -f $t.gz ]; then
			t=$t.gz
		fi	
		gunzip -cf $t | if [ -z "$(echo $r | pcut -t '|' -f 3)" ]; then
			cat
		else
			pcut -m 1 "_$(echo $r | pcut -t '|' -f 3)" 
		fi | rtool namerepl -r + | rtool flatten -s $HSTR | gzip -c > $DIR/compare/$(rname $r)/bregions.txt.gz
	fi

	if [ ! -f $DIR/compare/$(rname $r)/regions.txt.gz ]; then 
		# crop the foreground regions to background and include background regions as +
		echo "making $DIR/compare/$(rname $r)/regions.txt.gz"
		rtool flatten -s $HSTR $DIR/regions.txt.gz | rtool crop -s $HSTR - $DIR/compare/$(rname $r)/bregions.txt.gz | rtool merge - $DIR/compare/$(rname $r)/bregions.txt.gz | gzip -c > $DIR/compare/$(rname $r)/regions.txt.gz
	fi

	if [ ! -f $DIR/compare/$(rname $r)/regions_sizes.txt ]; then 
		# 1) region
		# 2) size of overlap with background
		# 3) number of regions
		# 4) number of regions that at least partially overlap with background
		# NOTE: for HSTR, regions that do not have a strand will cause weird things to happen. split these into separate + and - regions or exclude them!
		(
			grep-overlap -t$HSTR -p -s " " $DIR/regions.txt.gz $DIR/compare/$(rname $r)/bregions.txt.gz 
			gunzip -c $DIR/compare/$(rname $r)/regions.txt.gz | awk '{print "x", $0}'
		) | awk '
			$1 == "0" || $1 == "1"{Y[$2]++; Z[$2]+=$1}
			$1 == "x"{X[$2]+=$5-$4+1}
			
			END{
				Y["+"]
				for (x in Y)
					print x, X[x]+0, Y[x]+0, Z[x]+0
			}'  - > $DIR/compare/$(rname $r)/regions_sizes.txt 
	fi
done


########################### RUN MOTIF INSTANCE PIPELINE, IF NECESSARY 
if [ ! -z "$MOTIFFILE" ]; then
	OPTMMDIR=$DIR/optmm
	KMERSD=X

	# only use the first background region for motif matches
	r=${BACKREGIONS/ */}
	if [ ! -f $DIR/motifs.txt ]; then 
		if [ "$KMERS" != "Z" ]; then 
			# merge all the Kmers into one run
			rmc -N motifs-pval -v $DIR/motifs-pval.rmc -c "MotifAddPValCutoff.sh -k '$KMERS' -s '_%dmer' -m $KEEPWITHMAX" -i $MOTIFFILE -f 1 -e - -n $NUMNODES -o $DIR/motifs.txt
		else
			awk -vOFS="\t" '/^>/{$1 = $1 "_Zmer"}; 1' $MOTIFFILE > $DIR/motifs.txt
		fi
	fi

	if [ ! -d $OPTMMDIR/Xmer/${MATCHMASKS/ */}/matches ]; then 
		# when $LOWMODE=1, we use -F "" to quit after control motif creation
		# NOTE: -r used to point to the original background regions file, but now we use this to allow for filtering
		FindOptMMParams.sh -N $NUMNODES -M -- -m $DIR/motifs.txt -r $DIR/compare/$(rname $r)/bregions.txt.gz -d $OPTMMDIR/Xmer -S $HSTR -E $OPTMMDIR -z $CONFIDZ -h 3 -o $ORG -k "$MATCHMASKS" -f "$MATCHESNAMES" -F "$([ "$LOWMODE" = 1 ] && echo "" || echo "0 1")"
	fi
else
	echo "fill in values for variables if they are not already there"
	# fill in values for variables if they are not already there
	if [ -z "$KMERS" ]; then
		KMERS=$(find $OPTMMDIR -maxdepth 2 -name final-controls.txt | awk -F"/" '$(NF-1) ~ /mer$/{print substr($(NF-1),1,length($(NF-1))-3)}' | sort -un | awk '{L=L " " $0}; END{print substr(L,2)}')
	fi
	KMERSD=$KMERS
	if [ -z "$MATCHMASKS" ]; then
		MATCHMASKS=$(
			for K in $(echo $KMERSD); do 
				find $OPTMMDIR/${K}mer -maxdepth 2 -name matches | awk -F'/' '{print $(NF-1)}'
			done | sort -u | awk '{L=L " " $0}; END{print substr(L,2)}'
		)
	fi

	if [ -z "$MATCHESNAMES" ]; then
		MATCHESNAMES=$(
			for K in $(echo $KMERSD); do for M in $(echo $MATCHMASKS); do
				if [ -d $OPTMMDIR/${K}mer/$M/matches ]; then 
					find $OPTMMDIR/${K}mer/$M/matches -type f 
				fi
			done; done | awk -F"/" 'sub(/[.]gz$/, "", $NF){print $NF}' | sort -k1,1n -u | awk '{L=L " " $0}; END{print substr(L,2)}'
		)
	fi
fi

########################### COUNT OVERLAP BETWEEN MOTIFS AND REGIONS
# region_sizes, match-counts, match-region-counts (except possibly gzip) should be identical to what is produced by RunCompareAnalysis.sh
if [ "$LOWMODE" = 0 ]; then
	echo "COUNT OVERLAP BETWEEN MOTIFS AND REGIONS"
	for r in $(echo $BACKREGIONS); do  for K in $(echo $KMERSD); do for M in $(echo $MATCHMASKS); do  for n in $(echo $MATCHESNAMES); do
		if [ -f $OPTMMDIR/${K}mer/$M/matches/$n.gz ]; then 
			mkdir -p $DIR/compare/$(rname $r)/${K}mer/$M/$n
			if [ ! -f $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-counts.txt.gz ]; then 
				cat <<-EOF
				grep-overlap -t$HSTR -r1 $OPTMMDIR/${K}mer/$M/matches/$n.gz $DIR/compare/$(rname $r)/regions.txt.gz | pcut -t '|' -T ' ' -f 2.1 1.1 + | gzip -c > $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-counts.txt.gz

				EOF
			fi

			if [ ! -f $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-region-counts.txt.gz ]; then 
				if [ $PRODUCERCF = 1 ]; then 
					# NOTE: unique identifier in 6th column ensures we don't miss duplicates
					cat <<-EOF
					grep-overlap -t$HSTR -r1 -o $OPTMMDIR/${K}mer/$M/matches/$n.gz $DIR/compare/$(rname $r)/bregions.txt.gz | grep-overlap -t$HSTR -r2 $DIR/regions.txt.gz - | awk -F'|' '
						\$1 != Last{delete Seen; Last=\$1}
						{
							split(\$2, A, /[\t ]/)
							split(\$1, B, /[\t ]/)
							if (!(A[1] in Seen))
							{
								X[B[1] OFS A[1]]++
								Seen[A[1]]
							}
						}
						END{for (x in X) print x, X[x]}
						' | gzip -c > $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-region-counts.txt.gz

					EOF
				else
					echo -n | gzip -c > $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-region-counts.txt.gz
				fi
			fi	
		fi
	done; done; done; done | rmc -N compare -v $DIR/compare.rmc -n 0 -f 2 -e -
else
	echo "we match the motifs and compute counts simultaneously "
	# NOTE: $BACKREGIONS will always have one element, K will always be X, M will always be 0, and n will always be 0.0
	for r in $(echo $BACKREGIONS); do for K in $(echo $KMERSD); do for M in $(echo $MATCHMASKS); do for n in $(echo $MATCHESNAMES); do
		mkdir -p $DIR/compare/$(rname $r)/${K}mer/$M/$n

		# create matches and count them simultaneously... avoids creating matches files when we don't care about conservation
		# NOTE: names have to be mapped to _C##
		if [ ! -f $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-counts.txt.gz ]; then 
			rmc -N instcompare -v $DIR/compare.rmc -c "extract-mfa -k 1 $(GetOrgVars.sh $ORG ALIGNFILE) < %s | motif-match -n 1 -m $DIR/optmm/${K}mer/motifs-toscan.txt | grep-overlap -t$HSTR -r1 - $DIR/compare/$(rname $r)/regions.txt.gz  | pcut -t '|' -T ' ' -f 2.1 1.1 +" -i $DIR/optmm/regions.txt -e - -n $NUMNODES | pcut -f 1 2 +3 -T ' ' | awk '
				ARGIND==1{
					X[$1]=$1
					for (i=2; i<=NF; i++)
						X[$i] = $1 "_C" (i-1)
				}
				
				ARGIND==2{
					print $1, X[$2], $3
				}' $DIR/optmm/Xmer/final-controls.txt - | gzip -c > $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-counts.txt.gz
		fi

		# create empty match-region-counts.txt; this framework doesn't permit producing it
		if [ ! -f $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-region-counts.txt.gz ]; then 
			echo -n | gzip -c > $DIR/compare/$(rname $r)/${K}mer/$M/$n/match-region-counts.txt.gz
		fi
	done; done; done; done
fi

########################### PRODUCE ENRICHMENT OUTPUT FILE
# NOTE: 2-8 match ComputeMotifEnrich.awk
# 1)  Background / Kmer / MatchMask / Confidence
# 2)  Motif (stripped of _#kmer)
# 3)  Background region (both + and ++)
# 4)  Foreground region
# 5)  Background count of motif
# 6)  Foreground count of motif
# 7)  Background count of control motifs
# 8)  Foreground count of control motifs
# 9)  Background region size
# 10) Foreground region size
# 11) Number of regions in foreground
# 12) Number of regions with a motif match (*only when match-region-counts is produced)
# 13) $ZCORR log2 enrichment of motif relative to size of background
# 14) $ZCORR log2 enrichment of motif relative to controls
# 15) $ZCORR log2 odds relative to controls
# 16) -log10 p-value of enrichment relative to controls (negative for depleted)
# 17) Fraction of regions with motif|lower:upper (*only when match-region-counts is produced)
if [ ! -f $DIR/enrichments.txt.gz ]; then 
echo "working enrichment:",$DIR/enrichments.txt.gz
	for r in $(echo $BACKREGIONS); do for K in $(echo $KMERSD); do for M in $(echo $MATCHMASKS); do for n in $(echo $MATCHESNAMES); do
		fn=$DIR/compare/${r/|*/}/${K}mer/$M/$n
		if [ -f $fn/match-counts.txt.gz -a -f $fn/match-region-counts.txt.gz ]; then
			awk -vRF=$DIR/compare/${r/|*/}/regions_sizes.txt -vN="${r/|*/}/%s/$M/$n" -vBL="+ ++ $ENRICHBACKS" '
				{
					c = sub(/_C[0-9]*$/, "", $2)
					R[$1]
					M[$2]
					C[ARGIND,$1, $2, c] += $3
				}

				END{
					split(BL,A)
					for (i in A) if (A[i] != "") BR[A[i]]
					while (getline < RF) 
					{
						S[$1] = $2
						NS[$1] = $4
					}
					for (m in M) 
					{
						ms = gensub(/_[0-9Z]*mer$/, "", 1, m)
						kmer = (m ~ /_[0-9Z]*mer$/) ? gensub(/^.*_/, "", 1, m) : ""
						
						for (r in R) if (r != "+") for (br in BR)
							print sprintf(N, kmer), ms, br, r, int(C[1, br, m, 0]), int(C[1, r, m, 0]), int(C[1, br, m, 1]), int(C[1, r, m, 1]), int(S[br]), int(S[r]), int(NS[r]), int(C[2, r, m, 0])
					}
				}
			' <(gunzip -cf $fn/match-counts.txt.gz) <(gunzip -cf $fn/match-region-counts.txt.gz)
		fi
	done; done; done; done | estats -d %f -i ....ya..zb -o "* C($ZCORR,2)" | estats -d %f -i ....yazb -o "* C($ZCORR,2) L($ZCORR,2) A(,10)" | estats -d %f -i ..........ta -o "* f|f(-$ZCORR):f($ZCORR)" | gzip -c > $DIR/enrichments.txt.gz
fi

########################### MAKE PLOTS
MOTIFS=$(
	gunzip -cf $DIR/enrichments.txt.gz | awk -vML="$MOTIFLIST" '
	BEGIN{
		split(ML,A,/ /)
		for (i in A) MX[A[i]]
	}
	{split($2,A,/_/)}
	ML == "" || (A[1] in MX) || ((A[1] "_" A[2]) in MX){print $2}' | sort -u
)

for colt                               coltn         col coln            DM    in \
	"Enrichment vs. size of regions"   "(log2)"      13  raw_enrich      0        \
	"Enrichment vs. shuffle motifs"    "(log2)"      14  control_enrich  0        \
	"Log2 odds"                        ""            15  log_odds        1        \
	"Enrichment (depletion) p-value"   "-log10"      16  enrich_p_val    2        \
	"Region recovery"                  ""            17  region_recov    1        \

do
	if [ $DM = 0 ]; then
		DMC='sprintf("%.1f", 2^V[i])'
	elif [ $DM = 1 ]; then
		DMC='sprintf("%.1f", V[i])'
	elif [ $DM = 2 ]; then
		DMC='v = (V[i] < 0 ? "(" : "") "10%%Super%%" sprintf("-%d", int(abs(V[i]))) "%%/Super%%" (V[i] < 0 ? ")" : "")'
	fi

	for rbgt                                   rbg  rbgn           in \
		""                                     "+"  all_background    \
		"- restricted to considered regions"   "++" other_regions     \
		
	do 
		# one plot for each cutoff/enrichment procedure
		# show all motifs that are significant
		# exclude from analysis regions that are the respective background (to improve FDR)
		if echo $ONPLOTS | awk -vP=0 '{exit ! and($1,2^P)}'; then
			for r in $(echo $BACKREGIONS); do for K in $(echo $KMERS); do for M in $(echo $MATCHMASKS); do for n in $(echo $MATCHESNAMES); do
				N=${r/|*/}/${K}mer/$M/$n
				mkdir -p $DIR/plots/by-cutoff/$N
				gunzip -cf $DIR/enrichments.txt.gz | awk -vN=$N -vB=$rbg -vC=$col -vOFS="\t" -vSUBSEP="\t" '
					function abs(x) { return x < 0 ? -x : x }
					$1==N && $3==B && $4 != "++" && $9 != $10{
						print -abs($16) * log(10) + log(2), $2, $4, $C
					}' | FDR.sh -c 1 -s $'\t' -L 1 -a $FDRFILTER -n 0 | awk -vOFS="\t" -vSUBSEP="\t" -vNT="$colt%%Newline%%$rbgt" '
					function abs(x) { return x < 0 ? -x : x }
					{
						V[$2, $3] = $4
						if (abs($4) > M)
							M = abs($4)
					}
					END{
						print "", NT "|NameH"
						for (i in V)
							print i, (V[i] / (M > 0 ? M : 1)) "||" (V[i] == 0 ? "" : ('$DMC')) | "ArrangeTable.sh"
					}' | tee $DIR/plots/by-cutoff/$N/$coln-$rbgn.txt | Heatmap.awk -F"\t" | gzip -c > $DIR/plots/by-cutoff/$N/$coln-$rbgn.svgz
			done; done; done; done
		fi

		# one plot for each motif/enrichment procedure
		if echo $ONPLOTS | awk -vP=1 '{exit ! and($1,2^P)}'; then
			for m in $(echo $MOTIFS); do 
				N=${m//\//}

				mkdir -p $DIR/plots/by-motif/$N

				gunzip -cf $DIR/enrichments.txt.gz | awk -vMotif=$m -vKMERS=$KMERS -vMMS=$MATCHMASKS -vBGS=$BACKREGIONS -vMNS=$MATCHESNAMES -vB=$rbg -vC=$col -vOFS="\t" -vNT="$colt $rbgt" -vRG=$REGIONLIST '
					function abs(x) { return x < 0 ? -x : x }
					$3==B && $4 != "++" && $2 == Motif{
						V[$1, $4] = $C
						if (abs($C) > M)
							M = abs($C)

						split($1,A,/[/]/)
						HasMN[A[2],A[4]]
					}
					END{
						print "", NT "|NameH"
						nKMERS = split(KMERS, KMERSa, / /)
						nMMS = split(MMS, MMSa, / /)
						nBGS = split(BGS, BGSa, / /)
						nMNS = split(MNS, MNSa, / /)
						nRG = split(RG, RGa, / /)

						for (mm=1; mm<=nMMS; mm++) 
						{
							if (nMMS > 1)
								print "", MMSa[mm] "|NameB|20|"

							curk = 0
							for (kmer=1; kmer<=nKMERS; kmer++) if ((KMERSa[kmer]"mer", MNSa[1]) in HasMN)
							{
								if ((curk++))
									print "", "|Spacer||mm" mm "-" kmer

								print "", "4%%Super%%-"(KMERSa[kmer]) "%%/Super%%|NameB|40|mm" mm
								for (mn=1; mn<=nMNS; mn++) for (rg=1; rg<=nRG; rg++) for (bg=1; bg<=nBGS; bg++) 
									if ((KMERSa[kmer]"mer", MNSa[mn]) in HasMN)
									{
										i = BGSa[bg] "/" KMERSa[kmer] "mer/" MMSa[mm] "/" MNSa[mn] SUBSEP RGa[rg]
										print RGa[rg], MNSa[mn] "|||" KMERSa[kmer] "-" MMSa[mm], V[i] == "" ? "-" : ((V[i] / (M > 0 ? M : 1)) "||" (V[i] == 0 ? "" : ('$DMC')))
									}
								print "", "|Name|40|" mm "-" kmer
							}

							if (nMMS > 1)
							{
								print "", "|Name|20|" mm
								if (mm != nMMS)
									print "", "|Spacer||mm" mm
							}

						}

						if (nMMS > 1)
							print "", "Match mask|NameH|20|"
						print "", "PWM P-value|NameH|40|"
						print "", "Conservation%%Newline%%confidence|NameH|73|"
					}' | tee $DIR/plots/by-motif/$N/$coln-$rbgn.txt | Heatmap.awk -F"\t" -vMarginTop=20 -vMarginRight=100 | gzip -c > $DIR/plots/by-motif/$N/$coln-$rbgn.svgz
			done
		fi

		# cumulative plots when motifs match regions
		if echo $ONPLOTS | awk -vP=2 '{exit ! and($1,2^P)}'; then
			for r in $(echo $BACKREGIONS); do 
				N=${r/|*/}
				mkdir -p $DIR/plots/cumul/$N

				for tmode in merge notmerge; do 
					for mmode in mean max; do 
						gunzip -cf $DIR/enrichments.txt.gz | awk -vAF="$ALIASFILE" -vBG=${r/|*/} -vB=$rbg -vC=$col -vOFS="\t" -vMM=$mmode -vTM=$tmode '
							BEGIN{
								if (AF != "")
									while (getline < AF) 
									{
										split($1,A,/;/)
										for (i in A) 
											M[A[i]] = A[1]
									}
							}
						
							{
								split($1,A1,/[/]/)
								split($2,A2,/_/)
								split($4,A4,/_/)
							}

							AF == ""{M[A2[1]] = A2[1]; M[A4[1]] = A4[1]}

							A1[1] == BG && $3 == B && (A2[1] in M) && (A4[1] in M) && M[A2[1]] == M[A4[1]]{
								N = A1[2] "/" A1[3] "/" A1[4]
								Ns[N]
								Rs[$4]
								
								if (MM == "mean")
								{
									Val[N,$4] += $C
									NVal[N,$4]++
								}
								else if (!((N,$4) in Val) || (0+$C) > Val[N,$4])
								{
									Val[N,$4] = 0+$C
									NVal[N,$4] = 1
									Max[N,$4] = $2
								}
							}
							function tn(x) { return TM=="merge" ? M[gensub(/_.*$/, "", 1, x)] : x }
							END{

								for (n in Ns)
								{
									for (r in Rs) if ((n,r) in Val)
										nRs[tn(r)]++
									for (r in Rs) if ((n,r) in Val)	
										print n, r, Val[n,r]/NVal[n,r], length(nRs) * nRs[tn(r)], Max[n,r]
									delete nRs
								}
							}
							' | sort -k3,3n | awk -vOFS="\t" '{
								T[$1] += 1/$4
								print $3, $1, T[$1], $2, $5
							}' | sort -k2,2 -k1,1n | tee $DIR/plots/cumul/$N/$coln-$rbgn-combine_${mmode}-$tmode.txt | Plot.awk -vXTitle="$colt $coltn" -vYTitle="Cumulative fraction" -vML=" " -vYMin=0 -vYMax=1 | gzip -c > $DIR/plots/cumul/$N/$coln-$rbgn-combine_${mmode}-$tmode.svgz
					done
				done
			done
		fi

		if echo $ONPLOTS | awk -vP=3 '{exit ! and($1,2^P)}'; then
			mkdir -p $DIR/plots/dataset-motif-match

			gunzip -cf $DIR/enrichments.txt.gz | awk -vAF="$ALIASFILE" -vKMERS=$KMERS -vMMS=$MATCHMASKS -vBGS=$BACKREGIONS -vMNS=$MATCHESNAMES -vB=$rbg -vC=$col -vOFS="\t" -vNT="$colt $rbgt" -vRG=$REGIONLIST '
				BEGIN{
					if (AF != "")
						while (getline < AF) 
						{
							split($1,A,/;/)
							for (i in A) 
								M[A[i]] = A[1]
						}
				}
				
				{
					split($2,A2,/_/)
					split($4,A4,/_/)
				}

				AF == ""{M[A2[1]] = A2[1]; M[A4[1]] = A4[1]}

				function abs(x) { return x < 0 ? -x : x }

				$3==B && $4 != "++" && (A2[1] in M) && (A4[1] in M) && M[A2[1]] == M[A4[1]]{
					V[$1, $2, $4] = $C
					if (abs($C) > Mx)
						Mx = abs($C)

					if (!(($2,$4) in MListS))
					{
						MListS[$2,$4]
						MList[$4] = MList[$4] " " $2
					}

					split($1,A,/[/]/)
					HasMN[A[2],A[4]]

				}
				END{
					print "", NT "|NameH"
					nKMERS = split(KMERS, KMERSa, / /)
					nMMS = split(MMS, MMSa, / /)
					nBGS = split(BGS, BGSa, / /)
					nMNS = split(MNS, MNSa, / /)
					nRG = split(RG, RGa, / /)

					for (mm=1; mm<=nMMS; mm++) 
					{
						if (nMMS > 1)
							print "", MMSa[mm] "|NameB|20|"

						curk = 0
						for (kmer=1; kmer<=nKMERS; kmer++) if ((KMERSa[kmer]"mer", MNSa[1]) in HasMN)
						{
							if ((curk++))
								print "", "|Spacer||mm" mm "-" kmer

							print "", "4%%Super%%-"(KMERSa[kmer]) "%%/Super%%|NameB|40|mm" mm
							for (mn=1; mn<=nMNS; mn++) for (rg=1; rg<=nRG; rg++) for (bg=1; bg<=nBGS; bg++) 
							{
								nMotif = split(substr(MList[RGa[rg]],2), Motifs, / /)

								if (nMotif > 0)
								{
									print RGa[rg] "|NameBR|25"

									for (xx=1; xx<=nMotif; xx++)
									{
										i = BGSa[bg] "/" KMERSa[kmer] "mer/" MMSa[mm] "/" MNSa[mn] SUBSEP Motifs[xx] SUBSEP RGa[rg]

										print Motifs[xx] "|||" RGa[rg], MNSa[mn] "|||" KMERSa[kmer] "-" MMSa[mm], V[i] == "" ? "-" : ((V[i] / (Mx > 0 ? Mx : 1)) "||" (V[i] == 0 ? "" : ('$DMC')))
									}
									print "|NameR|25|" RGa[rg]
									print "|Spacer||" RGa[rg]
								}
							}
							print "", "|Name|40|" mm "-" kmer
						}

						if (nMMS > 1)
						{
							print "", "|Name|20|" mm
							if (mm != nMMS)
								print "", "|Spacer||mm" mm
						}

					}

					if (nMMS > 1)
						print "", "Match mask|NameH|20|"
					print "", "PWM P-value|NameH|40|"
					print "", "Conservation%%Newline%%confidence|NameH|73|"
				}' | tee $DIR/plots/dataset-motif-match/$coln-$rbgn.txt | Heatmap.awk -F"\t" -vMarginTop=20 -vMarginRight=100 | gzip -c > $DIR/plots/dataset-motif-match/$coln-$rbgn.svgz
		fi
	done
done

