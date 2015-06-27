.PHONY: clean

all: splitter barcodes barcodes2 analysis2

splitter: sff2fastq/sff.o splitter.o
	gcc -g -o splitter sff2fastq/sff.o splitter.o -lstdc++

splitter.o: splitter.cc toolbox.h
	gcc -g -I. -c splitter.cc

barcodes: barcodes.o
	g++ -g -o barcodes barcodes.o -lpthread

barcodes.o: barcodes.cc toolbox.h
	g++ -Wl,--no-as-needed -g -Wall -std=c++11 -c barcodes.cc

barcodes2: barcodes2.o
	g++ -g -o barcodes2 barcodes2.o -lpthread

barcodes2.o: barcodes2.cc toolbox.h
	g++ -Wl,--no-as-needed -g -Wall -std=c++11 -c barcodes2.cc

analysis2: analysis2.o
	g++ -g -o analysis2 analysis2.o -lpthread
clean:
	rm splitter barcodes barcodes2 analysis2 *.o
