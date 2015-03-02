# Sequence manipulating code

## Prerequisites:

1. [SFF file](http://dl.dropbox.com/u/68829208/test.sff.zip) + [adaptor list](http://dl.dropbox.com/u/68829208/ionXpress_barcode.txt).

2. [FASTQ reads](http://dl.dropbox.com/u/68829208/testReads.fastq.gz) for exploratory analysis

## Task 1:

The first task is performed by the [benchmark.sh](https://github.com/koskot77/sophia/blob/master/benchmark.sh) script.
You are welcome to download, examine, and run it on your machine.
The script prepares the working directory and runs the [splitter.cc](https://github.com/koskot77/sophia/blob/master/splitter.cc) program.

The [splitter.cc](https://github.com/koskot77/sophia/blob/master/splitter.cc) program reads the adaptor list (_readAdaptors_ function) 
and iterates over the records in the SFF file looking for matches with the adaptors. Then it opens a file for each of the encountered
adaptors and copies the related chunks of data from the input file. The _Nreads_ field from the main header is corrected for each of
the output files to acknowledge the number of records copied.

The workhorse algorithms are collected in the [toolbox.h](https://github.com/koskot77/sophia/blob/master/toolbox.h) file and serve
a relatively *fast* and *efficient* engine searching for patterns in the genetic sequences. This is achieved by constructing a
hash table (implemented within _LookUpTable_ class) that returns adaptor's index once it is given a sequence that matches any of
the adaptors from the list (this costs only few operations independent on the adaptor length). Any string comparisons are avoided
by converting every genetic sequence into an array of integers (implemented within _NumericSequence_ class). The conversion is 
outlined in the _sequence2number_ function.

The processing time for the whole task #1 is proportional to  _O_(Nreads) * _O_(readLength).

## Task 2:

The [analysis](https://www.dropbox.com/s/japvy2quk7kxnqk/task2.pdf) can be reproduced with the
[tast2.Rmd](https://github.com/koskot77/sophia/blob/master/task2.Rmd) R markdown document.

In the process an interactive _Shiny_ application was also developed [deployed](https://koskot77.shinyapps.io/sophia/).
The bottom plot is interactive (try pointing to the dots in the bottom plot) and allows to relate the pattern IDs with
the original sequences.
