#!/bin/bash

# move input file to enumeration on the cli, using ; for line change and : for same family (same line)

source settings.sh

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUTLIST="$2"
    shift
    ;;
    -b|--bamPathsFile)
    BAMPATHS="$2"
    shift
    ;;
    -c|--crossrefs)
    CROSSREFIDs="$2"
    shift
    ;;
    -e|--exonReads)
    EXONREADS="$2"
    shift
    ;;
    -p|--parametersLine)
    PARAMETERS="$2"
    shift
    ;;
    -o|--outputFile)
    OUTFILE="$2"
    shift
    ;;
    -s|--sets)
    SETP="$2"
    shift
    ;;
    -r|--refcutoff)
    XREFCUTOFF="$2"
    shift
    ;;

    *)
            # unknown option
    ;;
esac
shift
done


USAGEMSG="
Usage: filtering.sh -i <inputFilePath> [-b bamInputFile -e minSubjectsExonReads] [-c ID1;ID2;ID3] -p \"-r -c=0 -ref /data/woods/..\"
	-o outputFilePrefix


-i	Input file

	An input file should look like this:

	/data/woods/Annotated_VCFs/1.annot.tab	/data/woods/Annotated_VCFs/2.annot.tab
	/data/woods/Annotated_VCFs/3.annot.tab
	/data/woods/Annotated_VCFs/15.annot.tab	/data/woods/Annotated_VCFs/16.annot.tab

	i.	All individuals pertaining to the same family in a single line.
	ii.	Individual file names belonging to same family tab-separated in each line.
	iii.	Single individuals from different families should each be written in their own lines.

	The script will check that files exists.

-b	BAM file names input file (optional, requires -e specified.)

	This file needs to have the same format of the input, so that each annotation (annot.tab) file
	can be associated to a BAM file. So, the equivant to the example file shown above would be

	BAM1Path	BAM2Path
	BAM3Path
	BAM15Path	BAM16Path

	Where paths are complete file paths, as in the previous case. The script will check that files
	exists.

-c	Cross reference identifiers (optional)

	A semi colon separated list of identifiers for cross references. The identifiers are expanded to
	$ANNOTATED_VCFS_PATH/<id>.annot.tab, so you need to make sure that those files exist. The script 
	checks for existance, and will warn if they are missing. The script will add the -ref part to
	the parameters line.

-e	Exon Read (optional, required when using -b)

	Number of individuals of those available in the bam file paths where the mutation must occur at least.	

-p	Parameters line

	A string ENCLOSED IN DOUBLE QUOTES which is given to the runMutationFilter. See parameters
	options on Katie's documentation (page 16 on approximately). Missing double quotes can make
	this wrapper fail miserably.

-o	Output file name prefix, full path.

-s	Set of parameters

	Allows to choose for predefined set of parameters (overrides -p). Options are:
		\"1\" : Recessive inheritance pattern for consanguineous individuals (depth cutoff = 5).
		\"2\" : Dominant inheritance pattern. (requires option -c to be set).
		\"3\" : Recessive inheritance pattern. (requires -b option, but sets -e=2).
		\"4\" : Only use filters of inheritance pattern, mutation type, and SNP frequency.

-r	CrossRefs reference cut-off (optional, required when using -c )
"

if ! [ -z $SETP ]
then
	if [ "$SETP" == "1" ]
	then
		PARAMETERS="-r -c=0 -d=5"
	elif [ "$SETP" == "2" ] 	
	then
		PARAMETERS="-Ext –snpF=0.1"
		XREFCUTOFF=2
		if [ -z $CROSSREFIDs ]
		then
			echo "-s option 2 requires option -c for CrossRefs ids being set, exiting"
			exit 1
		fi
	elif [ "$SETP" == "3" ]
	then
		PARAMETERS="-r"
		EXONREADS=2
		if [ -z $BAMPATHS ]
		then
			echo "-s option 3 requires option -b for file with BAM paths being set, exiting"
			exit 1
		fi
	elif [ "$SETP" == "4" ]
	then
		PARAMETERS="-filters in;mt;snp –Mut_Type m"
	else
		echo "Option $SETP for -s not recognized, please see usage by invoking without arguments"
		exit 1
	fi
fi	


if [ -z $INPUTLIST ]
then
	echo "Missing input option (-i)"
	echo "$USAGEMSG" 
	exit 1
fi


if ! [ -e $INPUTLIST ]
then
	echo "Input file $INPUTLIST not found"
	echo "$USAGEMSG" 
	exit 1
fi

if [[ -z  $PARAMETERS ]]
then
	echo "Missing parameters option (-p)"
	echo "$USAGEMSG"
	exit 1
fi

if [ -z $OUTFILE ]
then
	echo "Missing out file prefix option (-o)"
	echo "$USAGEMSG"
	exit 1
fi


PRESEED=`sha256sum $INPUTLIST | awk '{ print $1 }' | base64 | rev | head -c5`
SEED=$PRESEED`date +%s`
GROUPID=`echo $SEED | sha256sum | base64 | rev | head -c10; echo`

mkdir -p $RUNJOBSPATH

RUNFILTEREXEC=$RUNJOBSPATH/runFiltering_$GROUPID.sh

echo "Starting script creation..."

echo "#!/bin/bash" > $RUNFILTEREXEC
echo "#PBS -N MutationFilter$GROUPID" >> $RUNFILTEREXEC
echo "#PBS -l walltime=8:00:00" >> $RUNFILTEREXEC

cat settings.sh >> $RUNFILTEREXEC

XREFDEST=$MUTFILTXREFsDIR/xref_$GROUPID

# Use CrossReference

# copy desired .annot.tab to Cross-references dir, and remove any previous files there.
# using the same directory will make separate instances to clash, so we need to set up a cross reference folder
# for each run. The cross reference path is given as an argument to runFilter shell script, so we can control it
# and produce a separate one per run.
if ! [ -z $CROSSREFIDs ]
then
	if [ -z $XREFCUTOFF ]
	then
		echo "Using parameter -c for CrossReferences requires -r for the cross reference cut-off to be set"
		echo "removing $RUNFILTEREXEC"
		rm $RUNFILTEREXEC
		echo "Exiting!"
		exit 1
	fi
	# Check that cross ref files exist in the $ANNOTATED_VCFS_PATH 
	xrefIDs=$(echo $CROSSREFIDs | tr ";" "\n")

	for xref in $xrefIDs
	do
    		if ! [ -e $ANNOTATED_VCFS_PATH/$xref.annot.tab ]
		then
			echo "Missing annotation file for cross reference: $ANNOTATED_VCFS_PATH/$xref.annot.tab"
			echo "removing $RUNFILTEREXEC"
			rm $RUNFILTEREXEC
			echo "Exiting!"
			exit 1
		fi
	done

	cat copyCrossReferences_TEMPLATE.sh | sed "s/CROSSREFSLOT/$CROSSREFIDs/g" | sed "s+XREFDESTSLOT+$XREFDEST+g" >> $RUNFILTEREXEC
	PARAMETERS=$PARAMETERS" -ref $XREFDEST -ref_cutoff=$XREFCUTOFF"
fi

#echo "Past cross refs part"

# Use location coverage check

# Check bams mentioned in file.
if ! [ -z $BAMPATHS ] && [ -e $BAMPATHS ]
then
	echo "Got intp BAM part"
	if [ -z $EXONREADS ]
	then
		echo "Option -e (Exon reads) needs to be specified (integer > 0) when providing BAM files"
		echo "removing $RUNFILTEREXEC"
		rm $RUNFILTEREXEC
		echo "Exiting!"
		exit 1
	fi

	echo "Past first BAM if"
	bamPaths2Check=$(cat $BAMPATHS | tr "\t" "\n")
	
	for bamPath in $bamPaths2Check;
	do
		if ! [ -e $bamPath ]
		then
			echo "Missing BAM file for coverage check mentioned in -b given file: $bamPath"
			echo "removing $RUNFILTEREXEC"
			rm $RUNFILTEREXEC
			echo "Exiting!"
			exit 1
		fi

		if ! [ -e $bamPath.bai ]
		then
			echo "Missing BAM bai index file: $bamPath.bai"
			echo "removing $RUNFILTEREXEC"	
			run $RUNFILTEREXEC
			echo "Exiting!"
			exit 1
		fi
	done

	PARAMETERS=$PARAMETERS" -bamRef $BAMPATHS -exonRead=$EXONREADS"

fi

#echo "Past BAM parts"

echo "java -cp $MUTATIONFILTERPATH/src/:$MUTATIONFILTERPATH/* com.example.mutationfilter.FilterProgram $INPUTLIST $OUTFILE \"$PARAMETERS\"" >> $RUNFILTEREXEC

echo "rm -rf $XREFDEST" >> $RUNFILTEREXEC

chmod u+x $RUNFILTEREXEC

echo "About to send job to the cluster..."
echo "For your records, this is job execution $GROUPID"
echo "The execution script for this job is: $RUNFILTEREXEC"
echo "This file will include all settings used."
echo "The result file should be found here: $OUTFILE"

qsub -q $LONGQUEUE -d $RUNJOBSPATH $RUNFILTEREXEC

echo "Past qsub part"
