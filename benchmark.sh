#!/bin/bash

echo "Installing prerequisites"

git clone https://github.com/koskot77/sophia.git
cd sophia

git clone https://github.com/indraniel/sff2fastq.git
wget http://dl.dropbox.com/u/68829208/ionXpress_barcode.txt
wget http://dl.dropbox.com/u/68829208/test.sff.zip
wget http://dl.dropbox.com/u/68829208/testReads.fastq.gz

unzip test.sff.zip
mv R_2012_01_27_18_17_55_user_PGM-11-HNPCC_long_2_Auto_PGM-11-HNPCC_long_2_11.sff input.sff

cd sff2fastq && make && cd ../ && make

./sff2fastq/sff2fastq -o test.fastq input.sff

echo "Running clustering:"
time ./barcodes2 -i testReads.fastq

echo "The most frequent barcode candidate is:"
./analysis 

sed -e 's|,| |g' output.csv | awk 'BEGIN{m=0;p=""} $1!~/ID/{if(m<$3){m=$3;p=$2}} END{print p" "m}'

