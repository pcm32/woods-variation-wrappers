#!/bin/bash

# TODO multiple hypothesis testing should be default
# TODO add exome variant server queries.

source settings.sh

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--identifiers)
    IDENTIFIERS="$2"
    shift 
    ;;   
    -b|--useBAMs)
    useBAMs=0
    ;;
    -c|--runCorrection) 
    RUNCORRECTION=1
    shift    
    ;;
    -f|--femalesCount)
    FEMALESCOUNT="$2"
    shift
    ;;
    -m|--malesCount)
    MALESCOUNT="$2"
    shift
    ;; 
    -o|--outputFile)
    OUTFILE="$2"
    shift
    ;;
    -x|--bed)
    UseBED=0
    ;;
    -p|--processors)
    PROCS="$2"
    shift
    ;;

    *)
            # unknown option
    ;;
esac
shift
done

USAGEMSG="
Usage: toughMums.sh -i identifiers [-b bamIdentifiers] -f femalesCount -m malesCount [-c] [-o <outputFilePath>] -p numProcessors --bed

-i	Identifiers, \"{id1,id2,{id3..idN}}\", where each of these are expected to be found in the $ANNOTATED_VCFS_PATH/id1.annot.tab
	This uses shell brace expansion to list/enumerate identifiers.
	Documentation for brace expansion can be found here http://wiki.bash-hackers.org/syntax/expansion/brace.

	For instance, for all ids from 1 to 100, but skipping numbers 54 and 56, you would issue

	\"{{1..53},55,{57..100}}\"

	If the identifiers need to be padded by, say, 3 zeros on the left and the prefix \"IND\", you would issue

	\"IND000{{1..53},55,{57..100}}\"

	For brace expansion, the identifiers string needs to be quoted or double quoted.

-b	(optional) Flag Use BAM files for same identifiers as provided in -i
	Files will be expected to be in $BAMSPATH/<identifier>.bam and $BAMSPATH/<identifier>.bam.bai
	The script will check that files exists.

	For this option, brace expansion (as explained for option -i) can be used as well.

-f	Females count : the number of females. A non negative integer.

-m	Males count   : the number of males. A non negative integer.

-c	Flag to activate the multiple hypothesis testing correction.

-o	(optional) Out file path. Defaults to $TOUGHMUMSRESULTS/<groupID>_toughmums_out.txt

-p	(optional) Number of processors to use on the cluster node. Defaults to 4. Used for BAM processing, so giving
	more processors than BAM files will have no effect beyond the number of BAM files. Too many concurrent 
	processes might produce too much disk stress, running slower. 

-x      (optional) Use bedtools flag, replaces BED/BAM handling implemented previously by BEDtools arithmetics (faster).
"

if [ -z $IDENTIFIERS ]
then
	echo "Missing identifiers option (-i)"
	echo "$USAGEMSG"
	exit 1
fi

if [ -z $FEMALESCOUNT ]
then
	echo "Missing females count option (-f)"
	echo "$USAGEMSG"
	exit 1
fi

if [ -z $MALESCOUNT ]
then
	echo "Missing males count option (-m)"
	echo "$USAGEMSG"
	exit 1
fi

re='^[0-9]+$'
if ! [[ $MALESCOUNT =~ $re ]] ; then
   echo "error: option -m for males count is not a number" 
   exit 1
fi

if ! [[ $FEMALESCOUNT =~ $re ]] ; then
   echo "error: option -f for females count is not a number" 
   exit 1
fi

source functions.sh

# calling this should set $GROUPID
computeGroupID $IDENTIFIERS


if [ -z $OUTFILE ]
then
	echo "Missing outfile definition"
	OUTFILE=$TOUGHMUMSRESULTS/$GROUPID"_toughmums_out.txt"
	echo "Output written to $OUTFILE"
fi

if [ -z $PROCS ]
then
	PROCS=4
fi

OUTFILE_BASENAME=`basename $OUTFILE`
OUTFILE_DIR=`dirname $OUTFILE`

TOUGHMUMSTEMP=$TOUGHMUMSRESULTS"/toughMums_"$GROUPID"_tmp"
mkdir -p $TOUGHMUMSTEMP

reBrace='\{'
if [[ $IDENTIFIERS =~ $reBrace ]] ; then
	IFS=" " read -a identsArray <<< `eval echo $IDENTIFIERS`
else
	IFS=';' read -a identsArray <<< "$IDENTIFIERS"
fi

TOTALINDIVIDUALS=${#identsArray[@]}
IDENTIFIERS_DEST=$TOUGHMUMSTEMP/identifiers_tm_$GROUPID.txt
touch $IDENTIFIERS_DEST

BAMIDENTS_DEST=""
if [ $useBAMs ]; then
	BAMIDENTS_DEST=$TOUGHMUMSTEMP/bam_identifiers_tm_$GROUPID.txt
	touch $BAMIDENTS_DEST
fi


for ident in "${identsArray[@]}"
do
	if ! [ -e $ANNOTATED_VCFS_PATH/$ident.annot.tab ]
	then
		echo "Missing annotation file for indentifier $ident : $ANNOTATED_VCFS_PATH/$ident.annot.tab not available"
		rm $IDENTIFIERS_DEST
		echo "Exiting!"
                exit 1
	fi
        echo "ID "$ident
	echo $ANNOTATED_VCFS_PATH/$ident.annot.tab >> $IDENTIFIERS_DEST

	if [ $useBAMs ]; then
		if ! [[ -e $BAMSPATH/$ident.bam  && -e $BAMSPATH/$ident.bam.bai ]]
		then
			echo "Missing BAM file or BAM index (.bam.bai) for indentifier $ident in $BAMSPATH"
			rm $BAMIDENTS_DEST
			rm $IDENTIFIERS_DEST
			exit 1
		fi

		echo $BAMSPATH/$ident.bam >> $BAMIDENTS_DEST

	fi

done

SUMFM=$(($MALESCOUNT+$FEMALESCOUNT))
#echo "Total individuals "$TOTALINDIVIDUALS
#echo "Sum total "$SUMFM

if [ "$TOTALINDIVIDUALS" -ne "$SUMFM" ]
then
	echo "Female and male count sum $(($MALESCOUNT+$FEMALESCOUNT)) is different from total individuals provided (number of ids) $TOTALINDIVIDUALS"
	rm $IDENTIFIERS_DEST
	exit 1
fi

QUEUE=$SHORTQUEUE

TOUGHMUMSEXEC=$RUNJOBSPATH"/toughMums_"$GROUPID".sh"
touch $TOUGHMUMSEXEC

if [ $useBAMs ]; then
	echo "#PBS -l walltime=08:00:00" >> $TOUGHMUMSEXEC
        echo "#PBS -l nodes=1:ppn=$PROCS" >> $TOUGHMUMSEXEC
fi

echo "source $WRAPPERDIR/settings.sh" >> $TOUGHMUMSEXEC
echo "perl $TOUGHMUMSPATH/getCohortCounts.pl $IDENTIFIERS_DEST > $TOUGHMUMSTEMP/cohortCounts.txt" >> $TOUGHMUMSEXEC
echo "echo \"Done getCohortCounts.pl\" 1>&2" >> $TOUGHMUMSEXEC
echo "cd $TOUGHMUMSPATH" >> $TOUGHMUMSEXEC
echo "bash $TOUGHMUMSPATH/compareAllChromosomes.sh $TOUGHMUMSTEMP/cohortCounts.txt $OUTFILE $FEMALESCOUNT $MALESCOUNT" >> $TOUGHMUMSEXEC
echo "echo \"Done compareAllChromosomes.sh\" 1>&2" >> $TOUGHMUMSEXEC


if [ $useBAMs ]; then
	echo "IFS=\$'\\n' read -d '' -r -a lines < $IDENTIFIERS_DEST" >> $TOUGHMUMSEXEC
	echo "IFS=\$'\\n' read -d '' -r -a bamLines < $BAMIDENTS_DEST" >> $TOUGHMUMSEXEC
	FORPART="
# Generate lists of locations to check in each individual
checkpoint=\`date +%s\`
parallel --gnu -P $PROCS '
	NAME=\$(basename {} .annot.tab)
	useBed=$useBED
	# we should include generateLocationsToCheck in the repo, and modify it to produce a bed file
	# of the locations already sorted (NAME.locs.txt)
	if [ \$useBed ]
	then
	  perl $TOUGHMUMSPATHSC/generateLocationsToCheck.pl {} $OUTFILE bedfile | sed \"s/^chrM\\(\\s\\)/chrMT\\1/\" | sed \"s/^chr\\S+gl/chrGL/\" | sed \"s/\\(^chrGL\\S+\\)_random/\\1/\" | sed \"s/\\(^chrGL\\S+\\)/\\1\\.1/\" | sort -u -k 1,1 -k2,2n > $TOUGHMUMSTEMP/\$NAME.locs.bed
  	else
	  perl $TOUGHMUMSPATH/generateLocationsToCheck.pl {} $OUTFILE > $TOUGHMUMSTEMP/\$NAME.locs.txt
	  sort $TOUGHMUMSTEMP/\$NAME.locs.txt | uniq > $TOUGHMUMSTEMP/\$NAME\"_sortRes\"
	  mv $TOUGHMUMSTEMP/\$NAME\"_sortRes\" $TOUGHMUMSTEMP/\$NAME.locs.txt
        fi
	' ::: \${lines[@]}
#done
echo \"Done generateLocationsToCheck.pl\" 1>&2
checkpoint2=\`date +%s\`
echo \"Seconds taken : \"\$((checkpoint2-checkpoint)) 1>&2

# check actual locations for sequencing
#for i in \"\${bamLines[@]}\"
#do
parallel --gnu -P $PROCS '
        NAME=\$(basename {} .bam)
	useBed=$useBED
	if [ \$useBed ]
	then
	  bamToBed -i {} | sed \"s/^chrM\\(\\s\\)/chrMT\\1/\" | sed \"s/^chr\\S+gl/chrGL/\" | sed \"s/\\(^chrGL\\S+\\)_random/\\1/\" | sed \"s/\\(^chrGL\\S+\\)/\\1\\.1/\" | sort -k1,1 -k2,2n | intersectBed -a $TOUGHMUMSTEMP/\$NAME.locs.bed -b stdin -v -sorted | awk '\\''{ print \$1\":\"(\$2+1) }'\\'' > $TOUGHMUMSTEMP/\$NAME.unsequenced.txt
          # we add one to the coordinate at the end to move from 0-based nucleotide index in bam/bed to 1-based nucleotide based used by the existing scripts 
	else
      	  python $MUTATIONFILTERPATH/ExonReads/reads_only_locs.py {} $TOUGHMUMSTEMP/\$NAME.locs.txt > $TOUGHMUMSTEMP/\$NAME.unsequenced.txt
        fi
	echo \"$TOUGHMUMSTEMP/\$NAME.unsequenced.txt\" >> $TOUGHMUMSTEMP/namesOfUnsequencedLocFiles.txt
	' ::: \${bamLines[@]}
#done
echo \"Done reads_only_locs\" 1>&2
checkpoint3=\`date +%s\`
echo \"Seconds taken : \"\$((checkpoint3-checkpoint2)) 1>&2

# update frequency counts 
# Original counts denote lower bounds on frequencies because assumes all locations not in file are those corresponding to ref alleles
perl $TOUGHMUMSPATH/updateCounts.pl $TOUGHMUMSTEMP/namesOfUnsequencedLocFiles.txt $TOUGHMUMSTEMP/cohortCounts.txt > $TOUGHMUMSTEMP/newCohortCounts.txt
echo \"Done updateCounts.pl\" 1>&2
cd $TOUGHMUMSPATH
bash $TOUGHMUMSPATH/compareAllChromosomes.sh $TOUGHMUMSTEMP/newCohortCounts.txt $OUTFILE_DIR/UpperBounds_$OUTFILE_BASENAME $FEMALESCOUNT $MALESCOUNT
echo \"Done compareAllChromosomes\" 1>&2
checkpoint4=\`date +%s\`
echo \"Seconds taken (update and compare): \"\$((checkpoint4-checkpoint3)) 1>&2
rm $TOUGHMUMSTEMP/namesOfUnsequencedLocFiles.txt
"

	echo "$FORPART" >> $TOUGHMUMSEXEC
fi

echo "Submitting job $GROUPID to cluster"
echo "Intermediate files can be found in $TOUGHMUMSTEMP"

chmod u+x $TOUGHMUMSEXEC
qsub -d $TOUGHMUMSTEMP -q $QUEUE $TOUGHMUMSEXEC
