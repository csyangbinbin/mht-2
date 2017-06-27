#!/bin/bash

#Constants
COMPILE_PATH="build"
RUN_PATH="/$COMPILE_PATH/src/mht"
TEST_PATH="/$COMPILE_PATH/test/unit_tests"

#INPUT OPTIONS
RUN_OPTIONS=$1
INPUT_NAME=$2
OUTPUT_NAME=$3

#If the COMPILE doesn't exist
if [ ! -d $COMPILE_PATH ]; then
	mkdir $COMPILE_PATH;
fi

#Make the file
cd $COMPILE_PATH;
cmake ../; make -j4 #2>> error_output.txt;
cd ../;

#Run the code
if [ $RUN_OPTIONS = 0 ]
then
	rm $OUTPUT_NAME;
	.$RUN_PATH $INPUT_NAME > zed.csv;
	cp $OUTPUT_NAME ~/devel/gmm_phd/data;
fi

##Run the code
if [ $RUN_OPTIONS = 1 ]
then
   while read -r line
   do
   	NAME="$line"

	rm $NAME.csv
	.$RUN_PATH $NAME > $NAME.csv;
	mv $NAME.csv $OUTPUT_NAME;

   done < "$INPUT_NAME"
fi

#Run the tests
if [ $RUN_OPTIONS = 2 ] 
then
	.$TEST_PATH;
fi
