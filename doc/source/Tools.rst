Tools
=====

stats.summary
^^^^^^^^^^^^^

.. code:: bash

    #!bash
    slamdunk stats.summary -o 201601_mRNA-SLAMseq_pulse-chase_mapping_summary.csv 
    -s snps/*.txt -m mapped/*_slamdunk_mapped.bam 
    -f filtered/*_slamdunk_mapped_filtered.bam -n ../samples.csv


stats.tcperreadpos
^^^^^^^^^^^^^^^^^^

.. code:: bash

    #!bash
    slamdunk stats.tcperreadpos -t 12 -o results/conversionperreadposition/ 
    -r /project/ngs/philipp/slamseq/ref/GRCm38.fa -l 55 -s snps/ filtered/*_slamdunk_mapped_filtered.bam

dump.reads
^^^^^^^^^^

Print all infos per read into a file.

.. code:: bash

    #!bash
    slamdunk dump.reads -t 12 -o results/readinfo/ -r /project/ngs/philipp/slamseq/ref/GRCm38.fa 
    -s snps/ngm/ filtered/ngm/*_slamdunk_mapped_filtered.bam

