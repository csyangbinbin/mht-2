#!/bin/bash

#Constants
COMPILE_PATH="build"
RUN_PATH="/$COMPILE_PATH/src/mht"
TEST_PATH="/$COMPILE_PATH/test/unit_tests"

#PDF OUTPUT
OUTPUT_NAME="kalman"
DOT_NAME="$OUTPUT_NAME.dot"
PDF_NAME="$OUTPUT_NAME.pdf"

#INPUT OPTIONS
RUN_OPTIONS=$1

#If the COMPILE doesn't exist
if [ ! -d $COMPILE_PATH ]; then
	mkdir $COMPILE_PATH;
fi

#Make the file
cd $COMPILE_PATH;
cmake ../; make -j4;
cd ../;

#Run the code
if [ $RUN_OPTIONS = 0 ]
then
	.$RUN_PATH;
	#dot -Tpdf $DOT_NAME -o $PDF_NAME;
	#evince $PDF_NAME &
fi

#Run the tests
if [ $RUN_OPTIONS = 1 ] 
then
	.$TEST_PATH;
fi
