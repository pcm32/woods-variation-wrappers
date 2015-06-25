function checkFileExistanceExit {
	fileToCheck=$1
	idOfFile=$2
	if ! [ -e $fileToCheck ]
	then
		echo "Could not find file for ID $idOfFile $fileToCheck"
		echo "Exiting"
		exit 1
	fi
}

function computeGroupID {
	IDENTIFIERS=$1
	PRESEED=`echo $IDENTIFIERS | sha256sum | awk '{ print $1 }' | base64 | rev | head -c5`
	SEED=$PRESEED`date +%s`
	GROUPID=`echo $SEED | sha256sum | base64 | rev | head -c10; echo`
}

