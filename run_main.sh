#!/bin/bash

#Constants
COMPILE_PATH="build"
RUN_PATH="/$COMPILE_PATH/src/mht"
OUTPUT_NAME="kalman"
DOT_NAME="$OUTPUT_NAME.dot"
PDF_NAME="$OUTPUT_NAME.pdf"

#If the COMPILE doesn't exist
if [ ! -d $COMPILE_PATH ]; then
	mkdir $COMPILE_PATH;
fi

#Make the file
cd $COMPILE_PATH;
cmake ../; make -j4;
cd ../;

#Run the code
.$RUN_PATH;

#Look at the cluster graph
dot -Tpdf $DOT_NAME -o $PDF_NAME;
evince $PDF_NAME &
