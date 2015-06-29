#!/bin/bash

# add ability to run with multiple files.
# use group identifiers

source settings.sh

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--identifiers)
    IDENTIFIERS="$2"
    shift 
    ;;   
    -g|--genesFile) 
    QUERYGENESFILE="$2"
    shift    
    ;;

    *)
            # unknown option
    ;;
esac
shift
done

USAGEMSG="
Usage: findMutation.sh -i identifiers -g genesFile

-i	Identifiers, \"{id1,id2,{id3..idN}}\", where each of these are expected to be found in the $ANNOTATED_VCFS_PATH/id1.annot.tab
	This uses shell brace expansion to list/enumerate identifiers.
	Documentation for brace expansion can be found here http://wiki.bash-hackers.org/syntax/expansion/brace.

	For instance, for all ids from 1 to 100, but skipping numbers 54 and 56, you would issue

	\"{{1..53},55,{57..100}}\"

-g	Query genes file. The file containing the gene names to 
"

if [ -z $IDENTIFIERS ]
then
	echo "Missing identifiers option (-i)"
	echo "$USAGEMSG"
	exit 1
fi

if [ -z $QUERYGENESFILE ]
then
	echo "Missing query gene files option (-g)"
	echo "$USAGEMSG"
	exit 1
fi

if ! [ -e $QUERYGENESFILE ]
then
	echo "File with query genes (-g) not found, please make sure that a full path is given"
	echo "Exiting!!"
	echo "$USAGEMSG" 
	exit 1
fi

source functions.sh

# this defines the $GROUPID
computeGroupID $IDENTIFIERS

reBrace='\{'
if [[ $IDENTIFIERS =~ $reBrace ]] ; then
	IFS=" " read -a identsArray <<< `eval echo $IDENTIFIERS`
else
	IFS=';' read -a identsArray <<< "$IDENTIFIERS"
fi

#$GENEFILTERRESPATH
MUTEXEC=$RUNJOBSPATH/runFindMutationsSession_$GROUPID.sh

echo "#!/bin/bash" > $MUTEXEC
echo "#PBS -d $GENEFILTERPATH" >> $MUTEXEC

echo "Results should be available at:"

for ident in "${identsArray[@]}"; do
	IDENTIFIER=$ident
	
	if ! [ -e $ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab ]
	then
		echo "Annotated VCF file not available in "$ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab
		echo "Exiting!!"
		echo "Usage: findMutations.sh <completePathToFileWithGenesForQuery> <Identifier>" 
		exit 1
	fi
	
	echo "perl findMutationsinGene.pl $ANNOTATED_VCFS_PATH/$IDENTIFIER.annot.tab $QUERYGENESFILE > $GENEFILTERRESPATH/$IDENTIFIER.mutations.txt" >> $MUTEXEC
	

	echo $GENEFILTERRESPATH/$IDENTIFIER.mutations.txt

done	

echo "rm $MUTEXEC" >> $MUTEXEC
echo "#" >> $MUTEXEC

echo ""
echo "Sending findMutationsGene.pl job to cluster:"
# cat $MUTEXEC
# be more explicit on results
chmod u+x $MUTEXEC
MUTEXECID=`qsub -q $SHORTQUEUE -j oe -o /dev/null $MUTEXEC`

