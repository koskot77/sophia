.PHONY: clean splitter

splitter: sff2fastq/sff.o splitter.o
	gcc -g -o test sff2fastq/sff.o splitter.o -lstdc++

splitter.o: splitter.cc toolbox.h
	gcc -g -I. -c splitter.cc

clean:
	rm test *.o
