#!/bin/bash

source settings.sh

USAGEMSG="
Usage: dataFormatting.sh \"{identifiers_in_brace_expasion}\"

	The identifiers of the sample is what completes the following path: $OLD_VCFS_PATH/<identifier>.snps.vcf and, optionally if available
	$OLD_VCFS_PATH/<identifier>.indels.vcf . If both files are present, they will be merged, and then variant effect predictor is run. 
	If only the first file is present, then the variant effect predictor is run only on this file. The script will detect whether 
	the VCF was created with a different version of the human assembly/annotation than hg19, and run liftOver if necessary.

	Note that the identifiers within braces are always double quoted. This part will respond to bash brace expansion:

	- For an enumeration from 1 to 10, \"{1..10}\"
	- For an enumeration from 1 to 10 but skipping 5, \"{{1..4},{6..10}}\"
	- Just for 1 and 7, \"{1,7}\"

	So, ranges are defined as {start..stop} and single element enumerations as {eX,eY,eZ} (note the comma and double point usage).

	All jobs are sent to the cluster. Log files with info an any errors will be left at:

	$FORMATTINGTEMP/<groupIdentifier>/<identifier>.data.format.log.{1,2,3}

	where <groupIdentifier> is given at the end of the execution of the script.

	If the path where the file <identifier>.snps.vcf is expected to be found, it can be changed in the setting.sh file, modifying var
		OLD_VCFS_PATH
	"

if [  $# -le 1 ]; then 
	echo "$USAGEMSG"
	exit 1
fi

reBrace='\{'
if [[ $1 =~ $reBrace ]] ; then
	IFS=" " read -a identsArray <<< `eval echo $1`
else
	IFS=';' read -a identsArray <<< "$1"
fi

PRESEED=`echo $1 | sha256sum | awk '{ print $1 }' | base64 | rev | head -c5`
SEED=$PRESEED`date +%s`
GROUPID=`echo $SEED | sha256sum | base64 | rev | head -c10; echo`
TEMPDIR=$FORMATTINGTEMP/$GROUPID
mkdir -p $TEMPDIR

for FILEPREFIX in "${identsArray[@]}"; do

	LOG=$TEMPDIR/$FILEPREFIX".dataFormat.log"

	if [ -z $FILEPREFIX ]
	then
		echo "Missing indentifier"
		echo "$USAGEMSG"
		exit 1
	fi
	
	if ! [ -e $OLD_VCFS_PATH/$FILEPREFIX.snps.vcf ]
	then
		echo "File $OLD_VCFS_PATH/$FILEPREFIX.snps.vcf doesn't exist, exiting, run without arguments for help."
		exit 1
	fi	
	
	WAITFOR=""
	CORRECTASSEMBLY=1
	
	echo "Working on "$FILEPREFIX
		
	# first see whether files need merging or not
	if [ -e $OLD_VCFS_PATH/$FILEPREFIX.snps.vcf ] && [ -e $OLD_VCFS_PATH/$FILEPREFIX.indels.vcf ]
	then
		echo "Separate snps and indels files found, merging them..."
		echo "Separate snps and indels files found, merging them..." > $LOG.1
		MERGEPROC=`qsub -q $SHORTQUEUE -d $OLD_VCFS_PATH -j oe -o $LOG.1 -v FILE_IN=$FILEPREFIX $FORMATTING_TOOLS_PATH/runMergeSnpsIndels.sh`
		# we should copy the output of the mergeVCFs.e to a log and delete them in the end.
		WAITFOR='-W depend=afterok:'$MERGEPROC
		grep '^##' $OLD_VCFS_PATH/$FILEPREFIX.snps.vcf | grep -q hg19
		CORRECTASSEMBLY=$?
	else
	        grep '^##' $OLD_VCFS_PATH/$FILEPREFIX.vcf | grep -q hg19
	        CORRECTASSEMBLY=$?	
	fi
	
	
	if ! [ "$CORRECTASSEMBLY" == "0" ]
	then
		echo "VCF created with different version to hg19, running lift over..."
		echo "VCF created with different version to hg19, running lift over..." >> $LOG.2
		LIFTPROC=`qsub -q $SHORTQUEUE $WAITFOR -d $OLD_VCFS_PATH -j oe -o $LOG.2 -v FILE_IN=$FILEPREFIX $FORMATTING_TOOLS_PATH/runLiftOver.sh`
		# add the liftOverJob files to a log and delete them
		WAITFOR='-W depend=afterok:'$LIFTPROC
	fi
	
	cd $VEP_DIR
	echo "Running variant effect predictor..."
	echo "Running variant effect predictor..." >> $LOG.3
	VEPPROC=`qsub -q $LONGQUEUE $WAITFOR -j oe -o $LOG.3 -v FILE_IN=$FILEPREFIX ./runRecipe.sh`
done

echo "For your reference the group identifier for this jon is $GROUPID"
echo "Temp directory with errors and information is $TEMPDIR"
echo "Results should end in $ANNOTATED_VCFS_PATH"


