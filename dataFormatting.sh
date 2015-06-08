#!/bin/bash

source settings.sh

USAGEMSG="
Usage: dataFormatting.sh identifier

	The identifier of the sample is what complete the following path: $OLD_VCFS_PATH/<identifier>.snps.vcf and, optionally if available
	$OLD_VCFS_PATH/<identifier>.indels.vcf . If both files are present, they will be merged, and then variant effect predictor is run. 
	If only the first file is present, then the variant effect predictor is run only on this file. The script will detect whether 
	the VCF was created with a different version of the human assembly/annotation than hg19, and run lift over if necessary.

	All jobs are sent to the cluster. Log files with info an any errors will be lest at:

	$OLD_VCFS_PATH/<identifier>.data.format.log.{1,2,3}

	if the path where the <identifier>.snps.vcf is expected to be found, it can be changed in the setting.sh file, modifying var
		OLD_VCFS_PATH
	"

FILEPREFIX=$1

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

LOG=$OLD_VCFS_PATH/$FILEPREFIX".dataFormat.log"

touch $LOG

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






