# SeqArc
SeqArc compresses next-generation sequencing data in [FASTQ](http://en.wikipedia.org/wiki/Fastq) format. 

As one of few tools capable of handling data of various omics, species, read length, sequence generation and other features, SeqArc strikes the best balance between compression ratio and performance. 

It is achieved by block-united processing. A fixed volume of reads will be considered as a block, and reads within a block are encoded altogether. Each read contains ID, base sequence and quality scores. The three parts are compressed independently, with the most suitable algorithms respectively.

ID is compressed with adaptive binning and range coding. 
Base sequence will be mapped to reference if index built, replaced with alignment information and encoded with range coder. The base sequence unable to map will be estimated with 16-order model and encoded with range coder. 
Quality scores are estimated with some pre-set pattern and encoded with range coder. For users interested in higher compression ratio and slight modification of Q-scores, the R-Block algorithms are equipped within SeqArc. 

Besides, the compression output of SeqArc is packed with encapsulation format, so different parts can be independently accessed. It also conveniently supports following iteration and replacement of algorithms.

## Usage

### Params

    To build HASH index:
        SeqArc -i <ref.fa>
    OR to build BWA index:
        SeqArc -q -i <ref.fa>

    To compress:
        SeqArc -c [options] [ref.fa] -1 <input_file> -2 [input_file2] -o <result.arc>
    To decompress:
        SeqArc -d [options] [ref.fa] <result.arc> -o [fastq_prefix]

    Options:
          -t INT       Thread num for multi-threading, default as [1]
          -q           Compress with BWA index (no need for decompression)
          -h           Get instructions of this software
          -s           Create and load ref index frome memory, will cause memory strain
          -l DOUBLE    Set lossy compression of quality scores, the recommended value is 1.15
          -I INT       Set max insert-size for PE alignment
          -f           Force overwrite target file
          -P INT       Pipe out decompression result
             1         Pipe out SE reads, or PE1 reads
             2         Pipe out PE2 reads
             3         Pipe out each pair of PE reads in order

          -1           Input file for compression, SE or PE1 
          -2           Input file for compression, PE2 
          -o           Output file name
          -p           Output to the directory of input

### Practice
&nbsp;

### I. Index Building

It is highly recommended to have the reference fasta corresponding to your targeted fastq file. If not, the compression can only be performed with entropy coding, which much lowers the compression performance.

The following command helps you to build a HASH index of fasta for alignment afterwards:

    SeqArc -i ref.fa

You can also build a BWA index (while its alignment speed is lower than HASH's) with:

    SeqArc -q -i ref.fa


### II. Compression

As index built, the easiest way to compress single-end file is:

    SeqArc -c ref.fa 1.fq test

And for paired-end files:

    SeqArc -c ref.fa 1.fq 2.fq test

In these 2 examples, a file named "test.arc" will be created to **current directory**.

Note: If the index was built with "-q", the compression parameters **must** contain "-q" as well.

If no index included, just delete the 'ref.fa' and compression will begin without sequence alignment.
&nbsp;

### III. Regular Decompression
No matter single-end or paired-end input, the command to decompress test.arc is:

    SeqArc -d ref.fa test.arc back

For SE, a file named "back.fastq" will be created to **current directory**.
For PE, 2 files named "back1.fastq" and "back2.fastq" will be created to **current directory**.

For compression with index, just delete the 'ref.fa'.

Note: If the "back" part is not given, the name(s) of decompression result will be **exactly the same as compression input**.
&nbsp;

### IV. Decompression to standard-output
If you want to directly pass the decompressed reads to downstream analysis tools, use the parameter 'p' and they will be delivered to stdout.

For single-end data:

    SeqArc -d -p 1 ref.fa test.arc

For paired-end data:

All PE1 reads will be output(no PE2 output):

    SeqArc -d -p 1 ref.fa test.arc

All PE2 reads will be output(no PE1 output):

    SeqArc -d -p 2 ref.fa test.arc

PE1 and PE2 reads will be output interleavedly:

    SeqArc -d -p 3 ref.fa test.arc

&nbsp;

### End
