Dunks
=====

map
^^^

.. code:: bash

    #!bash
    slamdunk map -r /project/ngs/philipp/slamseq/ref/GRCm38.fa -o mapped/ -5 4 -t 12 
    /project/ngs/philipp/slamdunk-analysis/veronika/raw-data/201601_mRNA-SLAMseq_pulse-chase/*.fq.gz

Maps all sample files to the reference genome using the SlamSeq default parameters for NextGenMap (-e -i 0.95 --slamseq 2) and converts them to sorted and indexed BAM files. Output files have the same name as the input files with the prefix "_slamdunk_mapped.bam".

For example::

    34507_An312_wt-2n_mRNA-slamseq-autoquant_24h-star-R1.fq.gz -> 
    34330_An312_wt-2n_mRNA-slamseq-autoquant_0h-R1.fq_slamdunk_mapped.bam


Parameters:

* -t is the number of threads to use.

* -5 is the number of bp that will be trimmed at the 5' end. 

* -o directory where output files will be written.


TODO:

* Switch to local alignments

* Add more NextGenMap parameters

filter
^^^^^^

.. code:: bash

    #!bash
    slamdunk filter -t 12 -o filtered/ -t 12 mapped/*_slamdunk_mapped.bam

Filters mapped BAM files by mapping quality. Output files have the same name as input files with the prefix "_filtered"
Input: mapped BAM files
Output: filtered BAM files

TODO:

* Filter by identity

* Filter by NM

* Filter by mapped length

snp
^^^

.. code:: bash

    #!bash
    slamdunk snp -t 12 -o snps/ -f /project/ngs/philipp/slamseq/ref/GRCm38.fa 
    filtered/*_slamdunk_mapped_filtered.bam

Calls SNPs from samples.
Input: mapped/filtered BAM files
Output: one vcf file per sample (same name as input with prefix "_snp")

Parameters:

* -c is Minimimum coverage to call variant (default: 10)

* -a Minimimum variant fraction variant (default: 0.8)

# T->C count per UTR

count
^^^^^

.. code:: bash

    #!bash
    slamdunk count -t 12 -o results/tcounts/ -r /project/ngs/philipp/slamseq/ref/GRCm38.fa 
    -b /project/ngs/philipp/slamseq/ref/mm10_ensembl_3UTRs_final.bed -l 55 
    -s snps/ filtered/*_slamdunk_mapped_filtered.bam


stats.rates
^^^^^^^^^^^

.. code:: bash

    #!bash
    slamdunk stats.rates -t 12 -o results/rates/ -r /project/ngs/philipp/slamseq/ref/GRCm38.fa 
    filtered/*_slamdunk_mapped_filtered.bam
