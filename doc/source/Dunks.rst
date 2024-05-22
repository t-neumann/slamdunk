Dunks
=====

map
---

The *map* dunk is used to map reads to a given genome using `NextGenMap's <http://cibiv.github.io/NextGenMap>`_ SLAMSeq alignment settings.

.. code:: bash

    slamdunk map [-h] -r <reference fasta> -o <output directory> [-5 <bp to trim from 5' end>]
                 [-n <Output up to N alignments per multimapper>] [-a <maximum number of 3' As before trimming] [-t <threads>]
                 [-q] [-e] [-i <sample index>] [-ss] files [files ...]
                
Input
^^^^^
===================  ==================================================================================================================
File                 Description
===================  ==================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**files**            Samplesheet (see :ref:`sample-file`) or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files (wildcard \* accepted).
===================  ==================================================================================================================

Output
^^^^^^
============================  ===========================================================================================================
File                          Description
============================  ===========================================================================================================
**Mapped BAM file**           One BAM file per sample containing the mapped reads. 
**Mapped bam index file**     The raw unmapped reads in fastq / BAM format (multiple read files can be run at once).
============================  ===========================================================================================================

Output files have the same name as the input files with the prefix "_slamdunk_mapped.bam".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq.gz -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam

Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-r**     x         The reference fasta file.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-5**               Number of bases that will be hard-clipped from the 5' end of each read.
**-a**               Maximum number of As at the 3' end of a read.
**-n**               The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads with multiple alignments of equal best scores).
**-t**               The number of threads to use for this dunk. NextGenMap runs multi-threaded, so it is recommended to use more threads than available samples.
**-q**               Deactivates NextGenMap's SLAMSeq alignment settings. Can be used to align plain QuantSeq data instead of SLAMSeq data.
**-e**               Switches to semi-global alignment instead of local alignment.
**-i**               Run analysis only for sample <i>. Use for distributin slamdunk analysis on a cluster (index is 0-based).
**-ss**              Output BAM while mapping. Slower, but uses less hard disk.  
**file**   x         Samplesheet (see :ref:`sample-file`) or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

filter
------

The *filter* dunk is used to filter the raw alignments from the *map* dunk using multiple quality criteria to obtain the final alignments for all subsequent analyses.

.. code:: bash

    slamdunk filter [-h] -o <output directory> [-b <bed file>] [-mq <MQ cutoff>] [-mi <identity cutoff>]
                    [-nm <NM cutoff>] [-t <threads>] bam [bam ...]
                    
Input
^^^^^
=======  =======================================================
File     Description
=======  =======================================================
**bam**  The raw mapped reads in BAM format from the *map* dunk.
**-b**   (optional) Bed file with coordinates of 3'UTRs.
=======  =======================================================

Output
^^^^^^
==============================  ===============================================================================================================
File                            Description
==============================  ===============================================================================================================
**Filtered BAM file**           One BAM file per sample containing the filtered reads.
**Filtered bam index file**     One BAM index file accompanying each filtered BAM file.
==============================  ===============================================================================================================

Output files have the same name as the input files with the prefix "_filtered".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam

Parameters
^^^^^^^^^^
=========  ========  =================================================================================================================================================================================
Parameter  Required  Description
=========  ========  =================================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-b**               BED-file containing coordinates for 3' UTRs.
**-mq**              Minimum mapping quality required to retain a read.
**-mi**              Minimum alignment identity required to retain a read.
**-nm**              Maximum number of mismatches allowed in a read.
**-t**               The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples.
**bam**    x         BAM file(s) containing the raw mapped reads (wildcard \* accepted).
=========  ========  =================================================================================================================================================================================

-------------------------------------------------------------------------------------------------------------------------------

snp
---

The *snp* dunk is used to call variants on the final filtered alignments of the *filter* dunk using `VarScan2 <http://dkoboldt.github.io/varscan/>`_. Any called T->C SNPs from this dunk will be excluded in the subsequent
analyses to reduce the false-positive number. 

.. code:: bash

    slamdunk snp [-h] -o <output directory> -r <reference fasta> [-c <coverage cutoff>]
                 [-f <variant fraction cutoff>] [-t <threads>] bam [bam ...]
                    
Input
^^^^^
=======  ==============================================================
File     Description
=======  ==============================================================
**bam**  The final filtered reads in BAM format from the *filter* dunk.
=======  ==============================================================

Output
^^^^^^
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
^^^^^^^^^^
=========  ========  ==================================================================================================================================================================
Parameter  Required  Description
=========  ========  ==================================================================================================================================================================
**-h**               Prints the help.
**-r**     x         The reference fasta file.
**-o**     x         The output directory where all output files of this dunk will be placed. 
**-c**               Minimum coverage to call a variant.
**-f**               Minimum variant fraction to call a variant.
**-t**               The number of threads to use for this dunk. VarScan2 runs multi-threaded, so it is recommended to use more threads than available samples.
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  ==================================================================================================================================================================

------------------------------------------------------

count
-----

The *count* dunk calculates all relevant numbers on statistics of SLAMSeq reads for each given 3' UTR. Central output will be *tcount* table
(see :ref:`tcount-file`).

.. code:: bash

     slamdunk count [-h] -o <output directory> [-s <SNP directory>] -r <reference fasta> -b <bed file> [-c <conversion threshold>]
                     [-m] [-l <maximum read length>] [-q <minimum base quality>] [-t <threads] bam [bam ...]
                     
**Note:** Since QuantSeq is a strand-specific assay, only sense reads will be considered for the final analysis!
                    
Input
^^^^^
=======  =============================================================================================
File     Description
=======  =============================================================================================
**bam**  The final filtered reads in BAM format from the *filter* dunk.
**-s**   (optional) The called variants from the *snp* dunk to filter false-positive T->C conversions.
**-b**   Bed file with coordinates of 3'UTRs.
=======  =============================================================================================

Output
^^^^^^
======================  ==============================================================================================================
File                    Description
======================  ==============================================================================================================
**Tcount file**         A tab-separated *tcount* file per sample containing the SLAMSeq statistics (see :ref:`tcount-file`).
**Bedgraph file**       A bedgraph file per sample showing the T->C conversion rate on each covered reference T nucleotide.
**cB file (optional)**  A comma-separated *cB* file per sample containing all of the T->C mutational information (see :ref:`cB-file`).
======================  ==============================================================================================================

Output files have the same name as the input files with the prefix "_tcount".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount.csv
  
Parameters
^^^^^^^^^^
=========  ========  ================================================================================================================================================================================
Parameter  Required  Description
=========  ========  ================================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-s**               The output directory of the *snp* dunk containing the called variants.
**-r**     x         The reference fasta file.
**-b**     x         BED-file containing coordinates for 3' UTRs.
**-l**               Maximum read length (will be automatically estimated if not set).
**-c**               Number of T->C conversions in a read required to count it as a "TC" read.
**-m**               Flag to additionally create a cB.csv file, compatible with mixture modeling.
**-q**               Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use for this dunk. This dunk runs single-threaded so the number of threads should be equal to the number of available samples.
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  ================================================================================================================================================================================

------------------------------------------------------

all
---

The *all* dunk is used to run an entire *slamdunk* run at once. It sequentially calls the *map*, *filter*, *snp* and *count* dunks and
provides parameters to keep full control over all dunks.

.. code:: bash

    slamdunk all [-h] -r <reference fasta> -b <bed file> [-fb <bed file>] -o <output directory> [-5 <bp to trim from 5' end>] [-a <maximum number of 3' As before trimming]
                 [-n <Output up to N alignments per multimapper>] [-t <threads>] [-q] [-e] [-m] [-mq <MQ cutoff>]
                 [-mi <identity cutoff>] [-nm <NM cutoff>] [-mc <coverage cutoff>] [-mv <variant fraction cutoff>] [-mts]
                 [-rl <maximum read length>] [-mbq <minimum base quality>] [-i <sample index>] [-ss] files [files ...]
                
Input
^^^^^
=================== ===================================================================================================================
File                 Description
=================== ===================================================================================================================
**Reference fasta** The reference sequence of the genome to map against in fasta format.
**-b**              Bed file with coordinates of 3'UTRs.
**file**            Samplesheet (see :ref:`sample-file`) or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files (wildcard \* accepted).
=================== ===================================================================================================================

Output
^^^^^^

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
^^^^^^^^^^
=========  ========  =====================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================
**-h**     x         Prints the help.
**-r**     x         The reference fasta file.
**-b**     x         BED-file containing coordinates for 3' UTRs.
**-fb**              BED-file used to filter multimappers.
**-o**     x         The output directory where all output files of this dunk will be placed.
**-5**               Number of bases that will be hard-clipped from the 5' end of each read **[map]**.
**-a**               Maximum number of A at the 3' end of a read.
**-n**               The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads with multiple alignments of equal best scores) **[map]**.
**-t**               The number of threads to use for this dunk. NextGenMap runs multi-threaded, so it is recommended to use more threads than available samples
**-q**               Deactivates NextGenMap's SLAMSeq alignment settings. Can be used to align plain QuantSeq data instead of SLAMSeq data **[map]**.
**-e**               Switches to semi-global alignment instead of local alignment **[map]**.
**-m**               Use 3'UTR annotation to filter multimappers **[filter]**.
**-mq**              Minimum mapping quality required to retain a read **[filter]**.
**-mi**              Minimum alignment identity required to retain a read **[filter]**.
**-nm**              Maximum number of mismatches allowed in a read **[filter]**.
**-mc**              Minimum coverage to call a variant **[snp]**.
**-mv**              Minimum variant fraction to call a variant **[snp]**.
**-cb**              Flag to additionally create a cB.csv file, compatible with mixture modeling.
**-mts**             Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted. **[count]**.
**-rl**              Maximum read length (will be automatically estimated if not set) **[count]**.
**-mbq**             Minimum base quality for T->C conversions to be counted **[count]**.
**files**  x         Samplesheet (see :ref:`sample-file`) or a list of all sample BAM/FASTA(gz)/FASTQ(gz) files (wildcard \* accepted).
=========  ========  =====================================================================================================================================================
