Usage
=====

*Slamdunk* is a modular analysis software. Modules are dubbed *dunks* and each dunk builds upon the results from the previous dunks.


.. .. image:: img/slamdunk-pipeline.png
..   :width: 800px

The idea is to always process all samples with one command line call of a given dunk. Output files typically get the same name as the input files with a certain prefix (e.g. "slamdunk_mapped").
The command line interface follows the "samtools/bwa" style. Meaning that all commands are available through the central executable/script *slamdunk* (located in the bin directory)

Available dunks are::

    'map', 'filter', 'snp', 'count', 'all'

Calling a module with --help shows all possible parameters text:

.. code:: bash

    slamdunk map --help
    
    usage: slamdunk map [-h] -r REFERENCEFILE -o OUTPUTDIR [-5 TRIM5] [-n TOPN]
                    [-t THREADS] [-q] [-l]
                    bam [bam ...]
   positional arguments:
         bam                   Bam file(s)

   optional arguments:
         -h, --help            show this help message and exit
         -r REFERENCEFILE, --reference REFERENCEFILE
                        Reference fasta file
         -o OUTPUTDIR, --outputDir OUTPUTDIR
                        Output directory for mapped BAM files.
         -5 TRIM5, --trim-5p TRIM5
                        Number of bp removed from 5' end of all reads.
         -n TOPN, --topn TOPN  Max. number of alignments to report per read
         -t THREADS, --threads THREADS
                        Thread number
         -q, --quantseq        Run plain Quantseq alignment without SLAM-seq scoring
         -l, --local           Use a local alignment algorithm for mapping.

The flow of *slamdunk* is to first map your reads, filter your alignments, call variants on your final alignments and use these to calculate conversion rates, counts and various
statistics for your 3'UTRs.

.. image:: img/slamdunk_flow.png
   :width: 600px

All steps create a log file that has the same name as the output file. Typically there is one log file per sample and task (makes parallel execution easier).
Command line output is limited to a minimum at the moment. If a sample is finished a "." is printed (very basic progress bar).