#!/usr/bin/python
import sys
import re

#perl -ane '$l++;if($l==1){$m++;}if($l==2){$head=substr $F[0],0,15; if($head=~/AAGAGGATTC/){print $head."\n";}}if($l==4){$l=0;}' input.fastq

def readBarCodes(filename):
    bc = []
    with open(filename,'r') as f:
        while True:
            pair = f.readline()

            if len(pair) == 0:
                break

            (name,code) = pair.split('\t')
            name = name.rstrip()
            code = code.rstrip()

            bc.append( code )

    return bc

barcodes = readBarCodes('ionXpress_barcode.txt')

count = {}

with open(sys.argv[1]) as fastq:
    while True:
        id1 = fastq.readline()
        seq = fastq.readline().rstrip()
        id2 = fastq.readline()
        qua = fastq.readline().rstrip()

        if len(seq) == 0:
            break

        for code in barcodes:
            if re.search(code, seq[0:15]):

                 if code not in count:
                      count[code] = 0

                 count[code] += 1

for code in count:
    print code," -> ",count[code]

