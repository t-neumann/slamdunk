Usage
=====

*Slamdunk* is a modular analysis software. Modules are dubbed *dunks* and each dunk builds upon the results from the previous dunks.


.. .. image:: img/slamdunk-pipeline.png
..   :width: 800px

The idea is to always process all samples with one command line call of a given dunk. Output files typically get the same name as the input files with a certain prefix (e.g. "slamdunk_mapped").
The command line interface follows the "samtools/bwa" style. Meaning that all commands are available through the central executable/script *slamdunk* (located in the bin directory)
Wildcard characters to supply multiple files are recognized throughout *slamdunk*.

Available dunks are::

    'map', 'filter', 'snp', 'count', 'all'

Calling a module with --help shows all possible parameters text:

.. code:: bash

   usage: slamdunk all [-h] -r REFERENCEFILE -b BED [-fb FILTERBED] -o OUTPUTDIR [-5 TRIM5]
                    [-a MAXPOLYA] [-n TOPN] [-t THREADS] [-q] [-e] [-m]
                    [-mq MQ] [-mi IDENTITY] [-nm NM] [-mc COV] [-mv VAR]
                    [-mts] [-rl MAXLENGTH] [-mbq MINQUAL] [-i SAMPLEINDEX]
                    [-ss]
                    files [files ...]

   positional arguments:
      files                 Single csv/tsv file (recommended) containing all
                            sample files and sample info or a list of all sample
                            BAM/FASTA(gz)/FASTQ(gz) files

   optional arguments:
      -h, --help            show this help message and exit
      -r REFERENCEFILE, --reference REFERENCEFILE
                            Reference fasta file
      -b BED, --bed BED     BED file with 3'UTR coordinates
      -fb FILTERBED, --filterbed FILTERBED
                            BED file with 3'UTR coordinates to filter multimappers
      -o OUTPUTDIR, --outputDir OUTPUTDIR
                            Output directory for slamdunk run.
      -5 TRIM5, --trim-5p TRIM5
                            Number of bp removed from 5' end of all reads
                            (default: 12)
      -a MAXPOLYA, --max-polya MAXPOLYA
                            Max number of As at the 3' end of a read (default: 4)
      -n TOPN, --topn TOPN  Max. number of alignments to report per read (default:
                            1)
      -t THREADS, --threads THREADS
                            Thread number (default: 1)
      -q, --quantseq        Run plain Quantseq alignment without SLAM-seq scoring
      -e, --endtoend        Use a end to end alignment algorithm for mapping.
      -m, --multimap        Use reference to resolve multimappers (requires -n >
                            1).
      -mq MQ, --min-mq MQ   Minimum mapping quality (default: 2)
      -mi IDENTITY, --min-identity IDENTITY
                            Minimum alignment identity (default: 0.95)
      -nm NM, --max-nm NM   Maximum NM for alignments (default: -1)
      -mc COV, --min-coverage COV
                            Minimimum coverage to call variant (default: 10)
      -mv VAR, --var-fraction VAR
                            Minimimum variant fraction to call variant (default:
                            0.8)
      -mts, --multiTCStringency
                            Multiple T>C conversion required for T>C read
      -rl MAXLENGTH, --max-read-length MAXLENGTH
                            Max read length in BAM file
      -mbq MINQUAL, --min-base-qual MINQUAL
                            Min base quality for T -> C conversions (default: 27)
      -i SAMPLEINDEX, --sample-index SAMPLEINDEX
                            Run analysis only for sample <i>. Use for distributing
                            slamdunk analysis on a cluster (index is 1-based).
      -ss, --skip-sam       Output BAM while mapping. Slower but, uses less hard
                            disk.
                            
The flow of *slamdunk* is to first map your reads, filter your alignments, call variants on your final alignments and use these to calculate conversion rates, counts and various
statistics for your 3'UTRs.

.. image:: img/slamdunk_flow.png
   :width: 600px

All steps create a log file that has the same name as the output file. Typically there is one log file per sample and task (makes parallel execution easier).
Command line output is limited to a minimum at the moment. If a sample is finished a "." is printed (very basic progress bar).