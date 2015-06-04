#!/bin/bash

source settings.sh

FILEPREFIX=$1
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






