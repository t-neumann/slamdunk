Usage
=====

The idea is to always process all samples with one command line call. Output files typically get the same name as the input files with a certain prefix (e.g. "slamdunk_mapped"). *Please, let us know what you think about this!*

.. image:: img/slamdunk-pipline.png

The command line interface follows the "samtools/bwa" style. Meaning that all commands are available through one central executable/script called slamdunk (located in the bin direcotry)

Available modules are::

    'map', 'filter', 'snp', 'dedup', 'count', 'stats.rates', 'stats.summary', 
    'stats.tcperreadpos', 'dump.reads', 'all' (not implemented yet)


Calling a module with --help shows all possible parameters text:

.. code:: bash

    bin/slamdunk snp --help
    usage: slamdunk snp [-h] -o OUTPUTDIR -f FASTA [-c COV] [-a VAR] [-t THREADS]
    bam [bam ...]

    positional arguments:
      bam                   Bam file(s)

    optional arguments:
      -h, --help            show this help message and exit
    -o OUTPUTDIR, --outputDir OUTPUTDIR
      Output directory for mapped BAM files.
    -f FASTA, --fasta FASTA
      Reference fasta file
    -c COV, --min-coverage COV
      Minimimum coverage to call variant
    -a VAR, --var-fraction VAR 
      Minimimum variant fraction variant
    -t THREADS, --threads THREADS
      Thread number



All steps create a log file that has the same name as the output file. Typically there is one log file per sample and task (makes parallel execution easier).
Command line output is limited to a minimum at the moment. If a sample is finished a "." is printed (very basic progress bar).
At the moment the python code is pretty slow. As soon as everything works as expected, we will try to optimise the most crucial tasks.

Input:
^^^^^^

Output:
^^^^^^^
