.PHONY: clean splitter

splitter: ../Sophia/sff2fastq/sff.genome.o splitter.o
	gcc -g -o test ../Sophia/sff2fastq/sff.genome.o splitter.o -lstdc++

splitter.o: splitter.cc toolbox.h
	gcc -D__GENOME__ -g -I. -c splitter.cc

clean:
	rm test *.o
