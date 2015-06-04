#!/bin/bash

source settings.sh

QUERYGENESFILE=$1
IDENTIFIER=$2

#cd $GENEFILTERPATH

if ! [ -e $ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab ]
then
	echo "Annotated VCF file not available in "$ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab
	echo "Exiting!!"
	echo "Usage: findMutations.sh <completePathToFileWithGenesForQuery> <Identifier>" 
	exit 1
fi

if ! [ -e $QUERYGENESFILE ]
then
	echo "File with query genes (second argument) not found, please make sure that a full path is given"
	echo "Exiting!!"
	echo "Usage: findMutations.sh <completePathToFileWithGenesForQuery> <Identifier>" 
	exit 1
fi

MUTEXEC=$GENEFILTERRESPATH/runFindMutationsFor_$IDENTIFIER.sh
echo "#!/bin/bash" > $MUTEXEC
echo "#PBS -d $GENEFILTERPATH" >> $MUTEXEC
echo "perl findMutationsinGene.pl $ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab $QUERYGENESFILE > $GENEFILTERRESPATH/$IDENTIFIER.mutations.txt" >> $MUTEXEC
echo "rm $MUTEXEC" >> $MUTEXEC
echo "#" >> $MUTEXEC

echo "Sending findMutationsGene.pl job to cluster:"
cat $MUTEXEC
chmod u+x $MUTEXEC
MUTEXECID=`qsub -q $SHORTQUEUE -j oe -o /dev/null $MUTEXEC`


