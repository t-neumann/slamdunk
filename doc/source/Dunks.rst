Dunks
=====

map
^^^

The *map* dunk is used to map reads to a given genome using `NextGenMap's <http://cibiv.github.io/NextGenMap>`_ SLAMSeq alignment settings.

.. code:: bash

    slamdunk map [-h] -r <reference fasta> -o <output directory> [-5 <bp to trim from 5' end>]
                 [-n <Output up to N alignments per multimapper>] [-t <threads>]
                 [-q] [-l] bam [bam ...]
                
Input
"""""

* **Reference fasta**: The reference sequence of the genome to map against in fasta format.


* **bam**: The raw unmapped reads in fastq / BAM format (multiple read files can be run at once).

Output
""""""

* **Mapped BAM file**: One BAM file per sample containing the mapped reads. Output files have the same name as the input files with the prefix "_slamdunk_mapped.bam".
   For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq.gz -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam
    
* **Mapped bam index file**: One BAM index file accompanying each mapped BAM file.

* **Mapped bam flagstat file**: One BAM flagstat file accompanying each mapped BAM file containing basic statistics on the mapped BAM file.

Parameters
""""""""""

* **-h** Prints the help.

* **-r** The reference fasta file.

* **-o** The output directory where all output files of this dunk will be placed. 

* **-5** Number of bases that will be hard-clipped from the 5' end of each read.

* **-n** The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads with multiple alignments of equal best scores).

* **-t** The number of threads to use for this dunk. NextGenMap runs multi-threaded, so it is recommended to use more threads than available samples.

* **-q** Deactivates NextGenMap's SLAMSeq alignment settings. Can be used to align plain QuantSeq data instead of SLAMSeq data.

* **-l** Switches to local alignment instead of semi-global alignment. Semi-global alignments are the default for NextGenMap.

* **bam** Fastq/BAM file(s) containing the raw unmapped reads. Can be multiple if multiple samples are analysed simultaneously.

filter
^^^^^^

The *filter* dunk is used to filter the raw alignments from the *map* dunk using multiple quality criteria to obtain the final alignments for all subsequent analyses.

.. code:: bash

    slamdunk filter [-h] -o <output directory> [-b <bed file>] [-mq <MQ cutoff>] [-mi <identity cutoff>]
                    [-nm <NM cutoff>] [-t <threads>] bam [bam ...]
                    
Input
"""""

* **bam**: The raw mapped reads in BAM format from the *map* dunk.

Output
""""""

* **Filtered BAM file**: One BAM file per sample containing the filtered reads. Output files have the same name as the input files with the prefix "_filtered".
   For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam
    
* **Filtered bam index file**: One BAM index file accompanying each filtered BAM file.

* **Filtered bam flagstat file**: One BAM flagstat file accompanying each filtered BAM file containing basic statistics on the mapped BAM file.


Parameters
""""""""""

* **-h** Prints the help.

* **-o** The output directory where all output files of this dunk will be placed. 

* **-b** BED-file containing coordinates for 3' UTRs.
     
     * This file will be used to filter multi-mapping reads:
     
         * Any reads having only alignments in one given 3'UTR and otherwise non 3'UTR regions, will be kept.
         * Any reads with alignments to multiple 3'UTRs will be discarded.

* **-mq** Minimum mapping quality required to retain a read.

* **-mi** Minimum alignment identity required to retain a read.

* **-nm** Maximum number of mismatches allowed in a read.

* **-t** The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples.

* **bam** Fastq/BAM file(s) containing the raw mapped reads. Can be multiple if multiple samples are analysed simultaneously.


snp
^^^

The *snp* dunk is used to call variants on the final filtered alignments of the *filter* dunk using `VarScan2 <http://dkoboldt.github.io/varscan/>`_. Any called T->C SNPs from this dunk will be excluded in the subsequent
analyses to reduce the false-positive number. 

.. code:: bash

    slamdunk snp [-h] -o <output directory> -f <reference fasta> [-c <coverage cutoff>]
                 [-a <variant fraction cutoff>] [-t <threads>] bam [bam ...]
                    
Input
"""""

* **bam**: The final filtered reads in BAM format from the *filter* dunk.

Output
""""""

* **VCF file**: One `VCF file <http://www.1000genomes.org/wiki/Analysis/vcf4.0/>`_ per sample containing the called variants. Output files have the same name as the input files with the prefix "_snp".
   For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_snp.vcf
  
Parameters
""""""""""

* **-h** Prints the help.

* **-f** The reference fasta file.

* **-o** The output directory where all output files of this dunk will be placed. 

* **-c** Minimum coverage to call a variant.

* **-a** Minimum variant fraction to call a variant.

* **-t** The number of threads to use for this dunk. VarScan2 runs multi-threaded, so it is recommended to use more threads than available samples.

* **bam** Fastq/BAM file(s) containing the final filtered reads. Can be multiple if multiple samples are analysed simultaneously.

count
^^^^^

The *count* dunk calculates all relevant numbers on statistics of SLAMSeq reads for each given 3' UTR. Central output will be *tcount* table.

.. code:: bash

     slamdunk count [-h] -o <output directory> [-s <SNP directory>] -r <reference fasta> -b <bed file> [-m]
                    -l <maximum read length> [-q <minimum base quality>] [-t <threads] bam [bam ...]
                    
Input
"""""

* **bam**: The final filtered reads in BAM format from the *filter* dunk.

* **-s**: (optional) The called variants from the *snp* dunk to filter false-positive T->C conversions.

Output
""""""

* **Tcount file**: A tab-separated *tcount* file per sample containing the SLAMSeq statistics. 
    
* **Bedgraph file**: A bedgraph file per sample showing the T->C conversion rate on each covered reference T nucleotide.

Output files have the same name as the input files with the prefix "_tcount".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount.csv
  
Parameters
""""""""""

* **-h** Prints the help.

* **-o** The output directory where all output files of this dunk will be placed. 

* **-s**: (optional) The output directory of the *snp* dunk containing the called variants.

* **-r** The reference fasta file.

* **-b** BED-file containing coordinates for 3' UTRs. For each entry in the BED-file the SLAMSeq statistics will be calculated.

* **-l** Maximum read length in the filtered BAM file.

* **-m** Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted.

* **-q** Minimum base quality for T->C conversions to be counted.

* **-t** The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples.

* **bam** Fastq/BAM file(s) containing the final filtered reads. Can be multiple if multiple samples are analysed simultaneously.

                    
