# SeqArc
SeqArc compresses next-generation sequencing data in [FASTQ](http://en.wikipedia.org/wiki/Fastq) format. Both entropy coding and alignment-based coding are adopted to achieve high compression ratio.

				
## Usage

### Params

	To build HASH index:
		SeqArc -i <ref.fa>
	OR to build BWA index:
		SeqArc -q -i <ref.fa>

	To compress:
		SeqArc -c [options] [ref.fa] <input_file> [input_file2] <compress_prefix>

	To decompress:
		SeqArc -d [options] [ref.fa] <***.arc> [fastq_prefix]
	   
		-t INT         Thread num for multi-threading, default as [1]
		-q             compress with BWA index (no need for decompression)
		-h             Get instructions of this software
		-s             create and load ref index frome memory, will cause memory strain
		-l double      Set Lossy Compression of Quality Scores, the best value is 1.15
		-I INT         Set PE Align Maxinsr
		-f             Force overwrite target file

### Compression
First of all, you need to have the reference fasta corresponding to your targeted fastq file. If not, the compression can only be performed with entropy coding, which much lowers the compression performance. And the command "SeqArc -i <ref.fa>" helps you to build a index of fasta for following alignment.

As index built, the easiest way to compress SE file is:

	SeqArc -c ref.fa 1.fq test
		
And for PE file:
	
	SeqArc -c ref.fa 1.fq 2.fq test

In these 2 conditions, a file named "test.arc" will be created.

### Decompression
No matter SE or PE input, the command to decompress \*.arc is:

	SeqArc -d ref.fa test back

For SE, a file named "basck.fastq" will be created.
For PE, 2 files named "back1.fastq" and "back2.fastq" will be created.
If the "back" part is not given, the name(s) of decompression result will be exactly the same as compression input. 
