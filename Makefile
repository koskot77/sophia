.PHONY: clean splitter barcode

all: splitter barcodes

splitter: sff2fastq/sff.o splitter.o
	gcc -g -o splitter sff2fastq/sff.o splitter.o -lstdc++

splitter.o: splitter.cc toolbox.h
	gcc -g -I. -c splitter.cc

barcodes: barcodes.o
	g++ -g -o barcodes barcodes.o -lpthread

barcodes.o: barcodes.cc toolbox.h
	g++ -Wl,--no-as-needed -g -Wall -std=c++0x -c barcodes.cc

clean:
	rm splitter barcode *.o
