Alleyoop
========

*Alleyoop* (**A**\ dditional s\ **L**\ amdunk he\ **L**\ p\ **E**\ r tools for an\ **Y** diagn\ **O**\ stics **O**\ r **P**\ lots) is a collection of tools for post-processing and running diagnostics on *slamdunk* analyses.
Similar to *slamdunk*, the command line interface follows the "samtools/bwa" style. Meaning that all commands are available through the central executable/script *alleyoop* (located in the bin directory).

**MultiQC:**

We implemented a module for `MultiQC <http://multiqc.info/>`_ to support summary reports of *Alleyoop* modules. Currently MultiQC supports the *summary*, *rates*, *utrrates*, *tcperreadpos* and *tcperutrpos* commands.

dedup
-----

This tool allows you to deduplicate a given bam-file. While many tools like `Picard tools <https://broadinstitute.github.io/picard/>`_ already collapses
reads with same start and end position on the chromosome, *alleyoop* only collapses reads with same start and end position, mapping to the same strand and identical
read sequence.

.. code:: bash

    alleyoop dedup [-h] -o <output directory> [-tc <# tc mutations>] [-t <threads>] bam [bam ...]
                
Input
^^^^^

===================  ==================================================================================================================================================================
File                 Description
===================  ==================================================================================================================================================================
**bam**              Fastq/BAM file(s) containing the final filtered reads from *slamdunk*. Can be multiple if multiple samples are analysed simultaneously (wildcard * is recognized).
===================  ==================================================================================================================================================================

Output
^^^^^^
============================  ===========================================================================================================
File                          Description
============================  ===========================================================================================================
**Deduplicated BAM file**     One BAM file per sample containing the deduplicated reads. 
============================  ===========================================================================================================

Output files have the same name as the input files with the prefix "_dedup".
For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered.bam -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_dedup.bam


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-tc**              Only select reads with -tc number of T>C mutations for deduplication.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.  
**bam**    x         Fastq/BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

collapse
--------

This tool allows you to collapse all 3'UTR entries of a *tcount* file into one single entry per 3'UTR (similar to the exons->gene relationship in a gtf-file).
All entries with identical 3'UTR IDs will be merged.

.. code:: bash

    alleyoop collapse [-h] -o <output directory> [-t <threads>] tcount [tcount ...]
                
Input
^^^^^

===================  ==================================================================================================================================================================
File                 Description
===================  ==================================================================================================================================================================
**Tcount file**      A tab-separated *tcount* file per sample containing the SLAMseq statistics (see :ref:`tcount-file`).
===================  ==================================================================================================================================================================

Output
^^^^^^
============================  =================================================================================================================================
File                          Description
============================  =================================================================================================================================
**Tcount file**               A tab-separated *tcount* file per sample containing the SLAMseq statistics (see :ref:`tcount-file`) with collapsed 3'UTR entries.
============================  =================================================================================================================================

Output files have the same name as the input files with the prefix "_collapsed".

For example::
   
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount.tsv -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped_filtered_tcount_collapsed.tsv


Parameters
^^^^^^^^^^
=============== ========  =====================================================================================================================================================================
Parameter       Required  Description
=============== ========  =====================================================================================================================================================================
**-h**                    Prints the help.
**-o**          x         The output directory where all output files of this tool will be placed.
**-t**                    The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.  
**Tcount file** x         A tab-separated *tcount* file per sample containing the SLAMseq statistics (see :ref:`tcount-file`) with collapsed 3'UTR entries.
=============== ========  =====================================================================================================================================================================

------------------------------------------------------  

rates
-----

This tool computes the overall conversion rates in your reads and plots them as a barplot.

.. code:: bash

    alleyoop rates [-h] -o <output directory> -r <reference fasta> [-mq <MQ cutoff>]
                   [-t <threads>] bam [bam ...]
                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**bam**              Fastq/BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**overallrates.csv**           Tab-separated table of the overall conversion rates. 
**overallrates.pdf**           Overall conversion rate plot file.
============================   ===========================================================================================================

Below is an example plot of the overall conversion rates of the reads in a sample. One can appreciate the typical excess of T->C conversion (A->G on minus strand)
of the SLAMseq technology for later labelling timepoints.

.. image:: img/stats.rates.png
   :width: 600px


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.  
**bam**    x         Fastq/BAM file(s) containing the final filtered reads. Can be multiple if multiple samples are analysed simultaneously (wildcard * is recognized).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

tccontext
---------

This tool computes the genomic context of all Ts in a read and plots them as barplot to inspect any biases in that direction.

.. code:: bash

    alleyoop tccontext [-h] -o <output directory> -r <reference fasta> [-mq <MQ cutoff>]
                       [-t <threads>] bam [bam ...]
                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**tccontext.csv**              Tab-separated table of the 5' and 3' T-contexts, separated by strand.
**tccontext.pdf**              T-context plot file.
============================   ===========================================================================================================

Below is an example plot of the T-context of all reads in a sample. On top you will find the 5' context of individual Ts, at the bottom the respective 3' context of the individual Ts.
Note, that these will not be reciprocal (see e.g. `this publication <http://www.sciencedirect.com/science/article/pii/S0888754305002600>`_).

.. image:: img/stats.TCcontext.png
   :width: 600px


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.  
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

utrrates
--------

This tool checks the individual conversion rates per 3'UTR and plots them as boxplots over the entire realm of 3'UTRs. Each conversion is normalized
to all possible conversions from it's starting base e.g. A->G / (A->A + A->G + A->C + A->T). 

.. code:: bash

    alleyoop utrrates [-h] -o <output directory> [-r <reference fasta>] [-mq <MQ cutoff>] [-m]
                      [-t <threads>] -b <bed file> -l <maximum read length> bam [bam ...]
                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-b**               Bed file with coordinates of 3'UTRs.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**mutationrates_utr.csv**      Tab-separated table with conversion reads, one UTR per line.
**mutationrates_utr.pdf**      UTR conversion rate plot file.
============================   ===========================================================================================================

Below is an example plot of conversion rates for all UTRs for a given sample. Typically, the individual conversions for a given starting
base are balanced and unbiased, except for T->C conversions in SLAMseq samples with longer labelling times. 

.. image:: img/stats.utrrates.png
   :width: 600px


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-m**               Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.
**-b**     x         Bed file with coordinates of 3'UTRs.
**-l**               Maximum read length in all samples (will be automatically estimated if not set).
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

snpeval
-------

This tool produces some QA about the quality of your variant calls: Ideally, your T>C SNP calls should be independently of the number
of reads with T>C conversions found in an UTR. Otherwise, this would mean that you call more T>C SNPs the more T>C reads you have and thus
you lose signal by falsely calling SNPs from true T>C conversions.

To assess this, the UTRs are ranked by the number of containing T>C reads and marked with a bar if also a T>C SNP was called in the respective UTR.
The list is then filtered for 3'UTRs with sufficient coverage to confidently call SNPs by using only the upper quartile of the 3'UTRs according to 
total read coverage.

The resulting plots will show once the distribution of SNPs across ranked 3'UTRs being blind to SNP information and including SNP information.
Ideally, one would see the SNPs biased towards 3'UTRs with high T>C read content in the blind situation and uniformly distributed across all 3'UTRs
in the informed situation.

This difference is also quantified using a GSEA-like Mann-Whitney-U test. 

.. code:: bash

    alleyoop snpeval [-h] -o <output directory> -s <SNP directory> -r <reference fasta> -b <bed file> [-c <coverage cutoff>]
                     [-f <variant fraction cutoff>] [-m] [-l <maximum read length>] [-q <minimum base quality>] [-t <threads>]
                     bam [bam ...]

Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-s**               Directory of called SNPs from *snp* dunk.
**-b**               Bed file with coordinates of 3'UTRs.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**SNPeval.csv**                Tab-separated table with read counts, T>C read counts and SNP indication, one UTR per line.
**SNPeval.pdf**                SNP evaluation plot.
============================   ===========================================================================================================

An example plot is coming soon!


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-s**     x         The output directory of the *snp* dunk containing the called variants.
**-r**     x         The reference fasta file.
**-q**               Minimum base quality for T->C conversions to be counted.
**-m**               Flag to activate the multiple T->C conversion stringency: Only T->C conversions in reads with more than 1 T->C conversion will be counted.
**-c**               Minimum coverage to call a variant.
**-f**               Minimum variant fraction to call a variant.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.
**-b**     x         Bed file with coordinates of 3'UTRs.
**-l**               Maximum read length in all samples (will be automatically estimated if not set).
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

summary
-------

This tool lists basic statistics of the mapping process in a text file.

.. code:: bash

    alleyoop summary [-h] -o <output file> [-t <directory of tcount files>] bam [bam ...]

Input
^^^^^

========================= =======================================================================================
File                      Description
========================= =======================================================================================
**bam**                   BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
**tcount file directory** (optional) Directory containing the associated tcount file(s) to the input BAM file(s).
========================= =======================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**outputfile**                 Tab-separated table with mapping statistics.
**outputfile_PCA.pdf**         PCA plot of the samples based on T>C read counts per UTR.
**outputfile_PCA.txt**         PCA values of the samples based on T>C read counts per UTR.
============================   ===========================================================================================================

The output file will be a tab-separated text file with the following columns:

============================   ===========================================================================================================
Column                         Content
============================   ===========================================================================================================
FileName                       Path to raw reads in BAM/fasta(gz)/fastq(gz) format.
SampleName                     Description of the sample.
SampleType                     Type of sample.  
SampleTime                     Timepoint of the sample in minutes.
Sequenced                      Number of sequenced reads.
Mapped                         Number of mapped reads.
Deduplicated                   Number of deduplicated reads.
Filtered                       Number of retained reads after filtering.
Counted                         Number of counted reads within UTRs **(optional: only if tcount file directory was supplied)**.
Annotation                     Annotation used for filtering.
============================   ===========================================================================================================

An example PCA plot is coming soon!

Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output file name.
**-t**               The directory of associated tcount file(s) to the supplied BAM file(s).
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

merge
-----

This tool merges *tcount* files of multiple samples into a single table based upon an expression of columns.

.. code:: bash

    alleyoop merge [-h] -o <output file> [-c <expression>] countFiles [countFiles ...]

Input
^^^^^

===================  =====================================================================================================
File                 Description
===================  =====================================================================================================
**countFiles**       A tab-separated *tcount* file per sample containing the SLAMseq statistics (see :ref:`tcount-file`).
===================  =====================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**outputfile**                 Tab-separated table merged *tcount* information based upon expression.
============================   ===========================================================================================================

Parameters
^^^^^^^^^^
============== ========  =====================================================================================================================================================================
Parameter      Required  Description
============== ========  =====================================================================================================================================================================
**-h**                   Prints the help.
**-o**         x         The output file name.
**-c**                   Column or expression used to summarize files (e.g. "TcReadCount / ReadCount")
**countFiles** x         A tab-separated *tcount* file per sample containing the SLAMseq statistics (see :ref:`tcount-file`).
============== ========  =====================================================================================================================================================================

------------------------------------------------------

tcperreadpos
------------

This tool calculates the individual mutation rates per position in a read treating T->C mutations separately. This plot can be used to search for biases
along reads. 

.. code:: bash

    alleyoop tcperreadpos [-h] -r <reference fasta> [-s <SNP directory>]
                          [-l <maximum read length>] -o <output directory> [-mq <MQ cutoff>]
                          [-t <threads>] bam [bam ...]
                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-s**               (optional) The called variants from the *snp* dunk to filter false-positive T->C conversions.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**tcperreadpos.csv**           Tab-separated table with mutation rates, one line per read position.
**tcperreadpos.pdf**           Plot of the mutation rates along the reads.
============================   ===========================================================================================================

Below is an example plot of mutation rates along all reads in a sample. Typically, one will see increasing error rates towards the end of a reads,
as for all Illumina reads. In addition, depending on how many bases were clipped from the 5' end of the reads, one will also observe higher error
rates in the beginning of the read as illustrated in the example plot. Finally, for SLAMseq samples with longer labelling times, the overall T->C 
conversions in the bottom plot will begin to increase compared to the overall background in the top plot.

.. image:: img/stats.tcperreadpos.png
   :width: 600px


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.
**-s**               The called variants from the *snp* dunk to filter false-positive T->C conversions.
**-l**               Maximum read length in all samples (will be automatically estimated if not set).
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

tcperutrpos
-----------

This tool calculates the individual mutation rates per position in an 3'UTR treating T->C mutations separately. This plot can be used to search for biases
along UTRs. Only most 3' 200bp of each UTR will be considered because: 
* Quantseq fragments are estimated have an average size of ~200bp
* This way, any dynamic binning biases can be avoided

.. code:: bash

   alleyoop tcperutrpos [-h] -r <reference fasta> -b <bed file> [-s <SNP directory>] 
                        [-l <maximum read length>] -o <output directory> [-mq <MQ cutoff>]
                        [-t <threads>] bam [bam ...]
                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-s**               (optional) The called variants from the *snp* dunk to filter false-positive T->C conversions.
**-b**               Bed file with coordinates of 3'UTRs.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**tcperutr.csv**               Tab-separated table with mutation rates, one line per UTR position.
**tcperutr.pdf**               Plot of the mutation rates along the UTRs.
============================   ===========================================================================================================

Below is an example plot of mutation rates along all UTRs in a sample. Typically, one will see increasing error rates towards the end of a UTRs.
For SLAMseq samples with longer labelling times, the overall T->C conversions in the bottom plot will begin to increase compared to the overall background in the top plot. 

.. image:: img/stats.tcperutrpos.png
   :width: 600px


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-b**     x         Bed file with coordinates of 3'UTRs.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.
**-s**               The called variants from the *snp* dunk to filter false-positive T->C conversions.
**-l**               Maximum read length in all samples (will be automatically estimated if not set).
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

------------------------------------------------------

dump
----

This tool outputs all available information calculated by *slamdunk* for each read in a sample.

.. code:: bash

   alleyoop dump [-h] -r <reference fasta> -s <SNP directory> -o <output directory>
                 [-mq <MQ cutoff>] [-t <threads>] bam [bam ...]

                
Input
^^^^^

===================  ===================================================================================================================================================================
File                 Description
===================  ===================================================================================================================================================================
**Reference fasta**  The reference sequence of the genome to map against in fasta format.
**-s**               The called variants from the *snp* dunk to filter false-positive T->C conversions.
**bam**              BAM file(s) containing the final filtered reads from *slamdunk* (wildcard \* accepted).
===================  ===================================================================================================================================================================

Output
^^^^^^
============================   ===========================================================================================================
File                           Description
============================   ===========================================================================================================
**readinfo.sdunk**             Tab-separated table with read info, one line per read
============================   ===========================================================================================================

The following columns are contained in the *readinfo* file:

============================   ===========================================================================================================
Column                         Description
============================   ===========================================================================================================
Name                           Name of the read
Direction                      Read was mapped on forward (1) or reverse (2) strand
Sequence                       Sequence of the read
Mismatches                     Number of mismatches in the read
tcCount                        Number of T->C conversion in the read
ConversionRates                List of all possible conversion in the read
============================   ===========================================================================================================


Parameters
^^^^^^^^^^
=========  ========  =====================================================================================================================================================================
Parameter  Required  Description
=========  ========  =====================================================================================================================================================================
**-h**               Prints the help.
**-o**     x         The output directory where all output files of this tool will be placed.
**-r**     x         The reference fasta file.
**-mq**              Minimum base quality for T->C conversions to be counted.
**-t**               The number of threads to use. All tools run single-threaded, so it is recommended to use as many threads as available samples.
**-s**     x         The called variants from the *snp* dunk to filter false-positive T->C conversions.
**bam**    x         BAM file(s) containing the final filtered reads (wildcard \* accepted).
=========  ========  =====================================================================================================================================================================

