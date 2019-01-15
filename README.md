# SeqArc
SeqArc compresses next-generation sequencing data in [FASTQ](http://en.wikipedia.org/wiki/Fastq) format. Both entropy coding and alignment-based coding are adopted to achieve high compression ratio.

##Installation
============
1. Clone the git repository,

		git clone https://github.com/BGI-flexlab/seqarc.git

2. Enter the src directory,
		
		cd src
			
3.	Compile,

		cmake . && make	
				
##Usage
=====
###Params


	To build index:
	  SeqArc -i <ref.fa>
	
	To compress:
	  SeqArc [options] [ref.fa] <input_file> [input_file2] <compress_prefix>
	    -l INT         min SMEM length to output [17]
	    -w INT         max interval size to find coordiantes [50]
	    -I INT         skip MEM mapped to over [-] places
	    -c INT         consider only the longest [2] sMEM for mapping
	    -E INT         drop out when getting a mapping result with [1] mismatch at most
	    -f INT         consider only the first [2] mapping result
	    -m INT         max mismatch to tolerate [3]
	    -s INT         max insert size between read1 and read2 [511]
	    -r INT         files with AlignRatio lower than [0.5] are processed with Fqzcomp only.
	
	    -S <level>     Sequence de novo compression level. 1-9 [3]
	                   Specifying '+' on the end (eg -s5+) will use
	                   models of multiple sizes for improved compression.
	    -N <level>     Quality compression level.  1-3 [2]
	    -n <level>     Name compression level.  1-2 [1]
	    -b             Use both strands in sequence hash table.
	    -e             Extra seq compression: 16-bit vs 8-bit counters.
	    -t INT         Thread num for multi-threading, default as [1]
	
	    -X             Enable generation/verification of check sums
	
	    -W             Show warning msg about abnormal base
	
	To decompress:
	   SeqArc -d [ref.fa] <compress_prefix> <fastq_prefix>

###Compression
First of all, you need to have the reference fasta corresponding to your targeted fastq file. If not, the compression can only be performed with entropy coding, which much lowers the compression performance. And the command "SeqArc -i \<ref.fa>" helps you to build a index of fasta for following alignment.

As index built, the easiest way to compress SE file is:

	SeqArc ref.fa 1.fq test
		
And for PE file:
	
	SeqArc ref.fa 1.fq 2.fq test

In these 2 conditions, a file named "test.arc" will be created.

###Decompression
No matter SE or PE input, the command to decompress *.arc is:

	SeqArc -d ref.fa test back

For SE, a file named "back.fastq"	 will be created.
For PE, 2 files named "back1.fastq" and "back2.fastq" will be created.


