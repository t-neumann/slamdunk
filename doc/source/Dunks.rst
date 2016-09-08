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
===================  ======================================================================================
File                 Description
===================  ======================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**bam**              The raw unmapped reads in fastq / BAM format (multiple read files can be run at once).
===================  ======================================================================================

Output
""""""
============================  ===========================================================================================================
File                          Description
============================  ===========================================================================================================
**Mapped BAM file**           One BAM file per sample containing the mapped reads. 
**Mapped bam index file**     The raw unmapped reads in fastq / BAM format (multiple read files can be run at once).
**Mapped bam flagstat file**  One BAM flagstat file accompanying each mapped BAM file containing basic statistics on the mapped BAM file.
============================  ===========================================================================================================

Output files have the same name as the input files with the prefix "_slamdunk_mapped.bam".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq.gz -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam

Parameters
""""""""""
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-r**     x         The reference fasta file.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-5**               Number of bases that will be hard-clipped from the 5' end of each read (default: 0).
**-n**               The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads with multiple alignments of equal best scores) (default: 1).
**-t**               The number of threads to use for this dunk. NextGenMap runs multi-threaded, so it is recommended to use more threads than available samples (default: 1).
**-q**               Deactivates NextGenMap's SLAMSeq alignment settings. Can be used to align plain QuantSeq data instead of SLAMSeq data.
**-l**               Switches to local alignment instead of semi-global alignment. Semi-global alignments are the default for NextGenMap.  
**bam**    x         Fastq/BAM file(s) containing the raw unmapped reads. Can be multiple if multiple samples are analysed simultaneously.
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

filter
^^^^^^

The *filter* dunk is used to filter the raw alignments from the *map* dunk using multiple quality criteria to obtain the final alignments for all subsequent analyses.

.. code:: bash

    slamdunk filter [-h] -o <output directory> [-b <bed file>] [-mq <MQ cutoff>] [-mi <identity cutoff>]
                    [-nm <NM cutoff>] [-t <threads>] bam [bam ...]
                    
Input
"""""
=======  =======================================================
File     Description
=======  =======================================================
**bam**  The raw mapped reads in BAM format from the *map* dunk.
**-b**   (optional) Bed file with coordinates of 3'UTRs.
=======  =======================================================

Output
""""""
==============================  ===============================================================================================================
File                            Description
==============================  ===============================================================================================================
**Filtered BAM file**           One BAM file per sample containing the filtered reads.
**Filtered bam index file**     One BAM index file accompanying each filtered BAM file.
**Filtered bam flagstat file**  One BAM flagstat file accompanying each filtered BAM file containing basic statistics on the filtered BAM file.
==============================  ===============================================================================================================

Output files have the same name as the input files with the prefix "_filtered".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam

Parameters
""""""""""
=========  ========  =================================================================================================================================================================================
Parameter  Required  Description
=========  ========  =================================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-b**               BED-file containing coordinates for 3' UTRs.
**-mq**              Minimum mapping quality required to retain a read (default: 2).
**-mi**              Minimum alignment identity required to retain a read (default: 0.8).
**-nm**              Maximum number of mismatches allowed in a read (default: -1).
**-t**               The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples (default: 1).
**bam**    x         BAM file(s) containing the raw mapped reads. Can be multiple if multiple samples are analysed simultaneously.
=========  ========  =================================================================================================================================================================================

-------------------------------------------------------------------------------------------------------------------------------

snp
^^^

The *snp* dunk is used to call variants on the final filtered alignments of the *filter* dunk using `VarScan2 <http://dkoboldt.github.io/varscan/>`_. Any called T->C SNPs from this dunk will be excluded in the subsequent
analyses to reduce the false-positive number. 

.. code:: bash

    slamdunk snp [-h] -o <output directory> -f <reference fasta> [-c <coverage cutoff>]
                 [-a <variant fraction cutoff>] [-t <threads>] bam [bam ...]
                    
Input
"""""
=======  ==============================================================
File     Description
=======  ==============================================================
**bam**  The final filtered reads in BAM format from the *filter* dunk.
=======  ==============================================================

Output
""""""
============  ===================================================================================================================
File          Description
============  ===================================================================================================================
**VCF file**  One `VCF file <http://www.1000genomes.org/wiki/Analysis/vcf4.0/>`_ per sample containing the called variants.
============  ===================================================================================================================

Output files have the same name as the input files with the prefix "_snp".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_snp.vcf
  
Parameters
""""""""""
=========  ========  ==================================================================================================================================================================
Parameter  Required  Description
=========  ========  ==================================================================================================================================================================
**-h**               Prints the help.
**-f**     x         The reference fasta file.
**-o**     x         The output directory where all output files of this dunk will be placed. 
**-c**               Minimum coverage to call a variant (default: 10).
**-a**               Minimum variant fraction to call a variant (default: 0.8).
**-t**               The number of threads to use for this dunk. VarScan2 runs multi-threaded, so it is recommended to use more threads than available samples (default: 1).
**bam**              BAM file(s) containing the final filtered reads. Can be multiple if multiple samples are analysed simultaneously.
=========  ========  ==================================================================================================================================================================

------------------------------------------------------

count
^^^^^

The *count* dunk calculates all relevant numbers on statistics of SLAMSeq reads for each given 3' UTR. Central output will be *tcount* table.

.. code:: bash

     slamdunk count [-h] -o <output directory> [-s <SNP directory>] -r <reference fasta> -b <bed file> [-m]
                     [-l <maximum read length>] [-q <minimum base quality>] [-t <threads] bam [bam ...]
                    
Input
"""""
=======  =============================================================================================
File     Description
=======  =============================================================================================
**bam**  The final filtered reads in BAM format from the *filter* dunk.
**-s**   (optional) The called variants from the *snp* dunk to filter false-positive T->C conversions.
**-b**   Bed file with coordinates of 3'UTRs.
=======  =============================================================================================

Output
""""""
==================  ===================================================================================================
File                Description
==================  ===================================================================================================
**Tcount file**     A tab-separated *tcount* file per sample containing the SLAMSeq statistics.
**Bedgraph file**   A bedgraph file per sample showing the T->C conversion rate on each covered reference T nucleotide.
==================  ===================================================================================================

Output files have the same name as the input files with the prefix "_tcount".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount.csv
  
Parameters
""""""""""
=========  ========  ================================================================================================================================================================================
Parameter  Required  Description
=========  ========  ================================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-s**               The output directory of the *snp* dunk containing the called variants.
**-r**     x         The reference fasta file.
**-b**     x         BED-file containing coordinates for 3' UTRs.
**-l**               Maximum read length (will be automatically estimated if not set).
**-m**               Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted.
**-q**               Minimum base quality for T->C conversions to be counted (default: 0).
**-t**               The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples (default: 1)
**bam**    x         BAM file(s) containing the final filtered reads. Can be multiple if multiple samples are analysed simultaneously.
=========  ========  ================================================================================================================================================================================

------------------------------------------------------

all
^^^

The *all* dunk is used to run an entire *slamdunk* run at once. It sequentially calls the *map*, *filter*, *snp* and *count* dunks and
provides parameters to keep full control over all dunks.

.. code:: bash

    slamdunk all [-h] -r <reference fasta> -b <bed file> -o <output directory> [-5 <bp to trim from 5' end>]
                 [-n <Output up to N alignments per multimapper>] [-t <threads>] [-q] [-l] [-m] [-mq <MQ cutoff>]
                 [-mi <identity cutoff>] [-nm <NM cutoff>] [-mc <coverage cutoff>] [-mv <variant fraction cutoff>] [-mts]
                 [-rl <maximum read length>] [-mbq <minimum base quality>] bam [bam ...]
                
Input
"""""
===================  ======================================================================================
File                 Description
===================  ======================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-b**               Bed file with coordinates of 3'UTRs.
**bam**              The raw unmapped reads in fastq / BAM format (multiple read files can be run at once).
===================  ======================================================================================

Output
""""""

One separate directory will be created for each dunk output:

==========  =============
Folder      Dunk
==========  =============
**map**     *map* 
**filter**  *filter* 
**snp**     *snp* 
**count**   *count* 
==========  =============

Parameters
""""""""""
=========  ========  =====================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================
**-h**     x         Prints the help.
**-r**     x         The reference fasta file.
**-b**     x         BED-file containing coordinates for 3' UTRs.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-5**               Number of bases that will be hard-clipped from the 5' end of each read (default: 0) **[map]**.
**-n**               The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads with multiple alignments of equal best scores) (default: 1) **[map]**.
**-t**               The number of threads to use for this dunk. NextGenMap runs multi-threaded, so it is recommended to use more threads than available samples (default: 1)
**-q**               Deactivates NextGenMap's SLAMSeq alignment settings. Can be used to align plain QuantSeq data instead of SLAMSeq data **[map]**.
**-l**               Switches to local alignment instead of semi-global alignment. Semi-global alignments are the default for NextGenMap **[map]**.
**-m**               Use 3'UTR annotation to filter multimappers **[filter]**.
**-mq**              Minimum mapping quality required to retain a read (default: 2) **[filter]**.
**-mi**              Minimum alignment identity required to retain a read (default: 0.8) **[filter]**.
**-nm**              Maximum number of mismatches allowed in a read (default: -1) **[filter]**.
**-mc**              Minimum coverage to call a variant (default: 10) **[snp]**.
**-mv**              Minimum variant fraction to call a variant (default: 0.8) **[snp]**.
**-mts**             Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted. **[count]**.
**-rl**              Maximum read length (will be automatically estimated if not set) **[count]**.
**-mbq**             Minimum base quality for T->C conversions to be counted (default: 0) **[count]**.
**bam**              Fastq/BAM file(s) containing the raw unmapped reads. Can be multiple if multiple samples are analysed simultaneously.
=========  ========  =====================================================================================================================================================

                    
