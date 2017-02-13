#!/bin/bash

# TODO find mutations that are present in the child (affected) but not in parents. Currently is producing the opposite result as required (showing currently the intersection of mutations between child and parents), when the result wanted are the new mutations shown by 
# the child that are not present in the parents. This should be available as a 4th option (-p) called de-novo dominant

source settings.sh

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUTLIST="$2"
    shift
    ;;
    -u|--unaffectedIDs)
    UNAFFECTEDIDS="$2"
    shift
    ;;
    -m|--motherID)
    MOTHERID="$2"
    shift
    ;;
    -d|--fatherID)
    FATHERID="$2"
    shift
    ;;
    -p|--inheritancePattern)
    INHERITANCEPAT="$2"
    shift
    ;;
    -o|--output-path)
    OUTPATH="$2"
    shift
    ;;
    -f|--familyIndex)
    FAMILYINDEX="$2"
    shift
    ;;

    *)
            # unknown option
    ;;
esac
shift
done


USAGEMSG="
Usage: furtherFiltering.sh -i <initialFilteringFileCompletePath> [-u ID1:ID2:ID3] [-m mother_identifier] [-d father_identifier] [-f family_index] -p inheritancePattern -o output_path


-i	Input file

	This should be the file produced by calling the filtering.sh routine. This is a full/complete path.

-u	Unaffected individuals Identifiers (optional)

	A colon separated list of identifiers for unaffected individuals. The identifiers are expanded to
	$ANNOTATED_VCFS_PATH/<id>.annot.tab, so you need to make sure that those files exist. The script 
	checks for existance, and will warn if they are missing. 

	For instance -u 1:4:7 would make the program expect the following files to be in place:

	$ANNOTATED_VCFS_PATH/1.annot.tab
	$ANNOTATED_VCFS_PATH/4.annot.tab
	$ANNOTATED_VCFS_PATH/7.annot.tab

-m	Mother identifier (optional)

	The identifier of the mother, expanded to $ANNOTATED_VCFS_PATH/<id>.annot.tab.

-d	Father (\"dad\") identifier (optional)

	The identifier of the father, expanded to $ANNOTATED_VCFS_PATH/<id>.annot.tab.

-f	Family index (within the input file), starts from zero for the first family in the input file.
	Only one family can be analyzed at a time. If non is given, the first family (zero) is used.

-p	Inheritance pattern

	The pattern encoded in one of the following options (one of them must be given):

	d   : dominant
	rnc : recessive, non-consanguineous
	rc  : recessive, consanguineous	

-o	Output path. Defaults to $CASECONTROLRESPATH.
"

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

MOTHERINPUT=""
if ! [ -z $MOTHERID ]
then
	if ! [ -e $ANNOTATED_VCFS_PATH/$MOTHERID.annot.tab ]
	then
		echo "Missing annotation file for mother: $ANNOTATED_VCFS_PATH/$MOTHERID.annot.tab"
		echo "Exiting!"
		exit 1
	fi

	MOTHERINPUT="--mother=$ANNOTATED_VCFS_PATH/$MOTHERID.annot.tab"
fi

FATHERINPUT=""
if ! [ -z $FATHERID ]
then
	if ! [ -e $ANNOTATED_VCFS_PATH/$FATHERID.annot.tab ]
	then
		echo "Missing annotation file for father: $ANNOTATED_VCFS_PATH/$FATHERID.annot.tab"
		echo "Exiting!"
		exit 1
	fi

	FATHERINPUT="--father=$ANNOTATED_VCFS_PATH/$FATHERID.annot.tab"
fi

if [ -z $FAMILYINDEX ]
then
	echo "Family index not set, using the first family (0)"
	FAMILYINDEX=0
fi
FAMILYINPUT="--family=$FAMILYINDEX"


INHERITANCE_INPUT=""
if ! [ -z $INHERITANCEPAT ]
then
	case $INHERITANCEPAT in
		d)
		INHERITANCE_INPUT="--dom"
		;;
		rnc)
		INHERITANCE_INPUT="--rec"
		;;
		rc)
		INHERITANCE_INPUT="--rec --con"
		;;

		*)
		echo "Inheritance pattern not recognized: $INHERITANCEPAT, see usage."
		echo "Exiting!"
		exit 1
	esac
else
	echo "Missing inheritance pattern"
	echo "$USAGEMSG"
	exit 1		
fi

if [ -z $OUTPATH ]
then
	OUTPATH=$CASECONTROLRESPATH
fi

PRESEED=`sha256sum $INPUTLIST | awk '{ print $1 }' | base64 | rev | head -c5`
SEED=$PRESEED`date +%s`
GROUPID=`echo $SEED | sha256sum | base64 | rev | head -c10; echo`

TEMP=$CASECONTROLTEMP/caseControl_$GROUPID
mkdir -p $TEMP

UNAFFDEST=$TEMP/unaffected.txt
touch $UNAFFDEST

# write unaffected exomes file for script
UNAFFECTED_INPUT=""
if ! [ -z $UNAFFECTEDIDS ]
then
	unafIDs=$(echo $UNAFFECTEDIDS | tr ":" "\n")

	for xref in $unafIDs
	do
    		if ! [ -e $ANNOTATED_VCFS_PATH/$xref.annot.tab ]
		then
			echo "Missing annotation file for unaffected individuals: $ANNOTATED_VCFS_PATH/$xref.annot.tab"
			echo "removing $UNAFFDEST"
			rm $UNAFFDEST
			echo "Exiting!"
			exit 1
		fi

		echo $ANNOTATED_VCFS_PATH/$xref.annot.tab >> $UNAFFDEST
	done
	UNAFFECTED_INPUT="--unAffected=$UNAFFDEST"
fi

RUNFFILTEREXEC=$RUNJOBSPATH/runFurtherFiltering_$GROUPID.sh
touch $RUNFFILTEREXEC

OUTFILENAME=$(basename $INPUTLIST)
OUTFILENAME=${OUTFILENAME%.*}

echo "perl $CASECONTROLPATH/caseControlCheck.pl --filterOut=$INPUTLIST $MOTHERINPUT $FATHERINPUT $UNAFFECTED_INPUT $INHERITANCE_INPUT $FAMILYINPUT > $OUTPATH/$OUTFILENAME-caseControl.txt" > $RUNFFILTEREXEC
echo "rm $UNAFFDEST" >> $RUNFFILTEREXEC

echo "Submitting further filtering job $GROUPID to the cluster"
cat $RUNFFILTEREXEC
echo "Temporary and log files will be in $TEMP"
echo "Results can be found in $OUTPATH/$OUTFILENAME*"

chmod u+x $RUNFFILTEREXEC

qsub -q short -d $TEMP $RUNFFILTEREXEC


