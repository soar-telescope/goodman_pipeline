#!/bin/bash

SOURCE_DIR=$(pwd)/goodman_pipeline/data/dcr_source/dcr

BINARY_FOLDER=$VIRTUAL_ENV/bin

if [ -d "$BINARY_FOLDER" ]; then
  echo "Installing dcr binaries in: $BINARY_FOLDER"
  if [ -d "$SOURCE_DIR" ]
	then
		echo "Compiling dcr"
		make --directory "$SOURCE_DIR"
		if [ -f "$SOURCE_DIR"/dcr ]
		then
			echo "Adding executing permission"
			chmod +x "$SOURCE_DIR"/dcr

#			if [ -d "$BINARY_FOLDER" ] &&  echo "$PATH" | grep -q "$BINARY_FOLDER"
#
#			then
#				echo "Copying binary file"
#				cp -v "$SOURCE_DIR"/dcr "$BINARY_FOLDER"
#			else
#				echo "Directory does not exist: $ENV_PATH"
#			fi
#
#			if [ -f "$BINARY_FOLDER"/dcr ]
#			then
#			    echo "DCR Installed Successfully on:  $BINARY_FOLDER"
#			else
#			    echo "DCR installation failed"
#			fi
		else
			echo "Unable to find compiled file"
		fi
	else
	    echo "Directory $SOURCE_DIR Does not exist"
	fi
else
  echo "Binary folder ${BINARY_FOLDER} does not exist"
fi
