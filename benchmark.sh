#!/bin/bash

echo "Installing prerequisites"

git clone https://github.com/indraniel/sff2fastq
wget http://dl.dropbox.com/u/68829208/ionXpress_barcode.txt
wget http://dl.dropbox.com/u/68829208/test.sff.zip
unzip test.sff.zip
mv R_2012_01_27_18_17_55_user_PGM-11-HNPCC_long_2_Auto_PGM-11-HNPCC_long_2_11.sff input.sff

cd sff2fastq && make && cd ../

./sff2fastq/sff2fastq -o data.txt input.sff

echo "Running a simplest pattern search:"
time for i in `cat ionXpress_barcode.txt | awk '{print $2}'`; do echo $i; grep $i data.txt | wc; done

echo "Running my program:"
time ./test

