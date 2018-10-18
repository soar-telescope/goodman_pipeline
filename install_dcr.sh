#!/bin/bash

SOURCE_DIR=$(pwd)/goodman_pipeline/data/dcr-source/dcr



if conda info | grep 'active environment' | grep -q 'None'

then
	echo 'You do not have any virtual environment activated'

else
	ENV_NAME=$(conda info | grep 'active environment' | sed 's/\<active environment\>//g' | sed "s/[: ]//g")
	ENV_PATH=$(conda info | grep 'active env location' | sed 's/\<active env location\>//g' | sed "s/[: ]//g")
	echo "Using Virtual Environment: " $ENV_NAME
	if [ -d $SOURCE_DIR ]
	then
		echo "Compiling dcr"
		make --directory $SOURCE_DIR
		if [ -f $SOURCE_DIR/dcr ]
		then
			echo "Adding executing permission"
			echo $SOURCE_DIR
			echo $PATH
			chmod +x $SOURCE_DIR/dcr

			if [ -d $ENV_PATH/bin ] &&  echo $PATH | grep -q $ENV_PATH/bin 
			
			then
				echo "Copying binary file"
				cp -v $SOURCE_DIR/dcr $ENV_PATH/bin
			else
				echo "Directory does not exist: " $ENV_PATH
			fi
		else
			echo "Unable to find compiled file"
		fi
	else
	    echo "Directory" $SOURCE_DIR "Does not exist"
	fi
fi

