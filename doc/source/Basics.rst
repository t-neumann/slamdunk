Basics
======

.. _conversionRates:

Conversion rates
----------------

The individual measured conversion rates within a 3' UTR are highly dependent on the T-content of a given UTR and the number of reads that cover those Ts.
The more reads containing Ts, the more accurate the overall estimation of the conversion rate of T>Cs in this given 3' UTR will be.

We have adressed this by using a 3' UTR positions based view instead of a read based view: For each T-position in a given 3' UTR, we record
the number of reads covering this position and the number of reads with T>C conversion at this position - the more reads, the more accurate the T>C conversion rate estimation will be.

Both values are summed over all T-positions in the UTR - this way, less reliable T-positions with little reads will have lower impact on the overall T>C conversion rate estimation than
T-positions with many reads.   


.. image:: img/conversionrates.png
   :width: 800px
   
**Note:** Since QuantSeq is a strand-specific assay, only sense reads will be considered for the final analysis!

Multimapper reconciliation
--------------------------

Definition and importance
^^^^^^^^^^^^^^^^^^^^^^^^^

Multimappers are reads which align to more than one position in the genome with equal best alignment scores.

Such reads make up a substantial part of the alignment and are generally ignored, but can provide valuable information if they can be unequivocally
reassigned to a given location, enhancing the overall signal and thus heavily influencing the overall accuracy.

However, there are also several caveats: For genes with isoforms or homologous genes, simply taking all alignments into account will heavily distort
the signal, since the reads will likely not account equally for all regions, but rather originate mainly from a single region which cannot be unambiguously assigned. 

Reassignment strategy
^^^^^^^^^^^^^^^^^^^^^

We have settled for a very conservative multimapper reassignment strategy:

Since the QuantSeq technology specifically enriches for 3' UTRs, we only consider alignments to annotated 3' UTRs supplied to *slamdunk* as relevant.

Therefore, any multimappers with alignments to a single 3' UTR and non-3'UTR regions (i.e. not annotated in the supplied reference) will be unequivocally
assigned to the single 3'UTR. If there are multiple alignments to this single 3'UTR, one will be chosen at random.

For all other cases, were a read maps to several 3' UTRs, we are unable to reassign the read uniquely to a given 3'UTR and thus discard it from the analysis.

In short, the procedure is as follows, illustrated by the example below:

#. Create an `Interval Tree <https://pypi.python.org/pypi/intervaltree/2.0.4>`_ from the supplied 3' UTR annotation
#. Check the overlaps of all multimapping reads with 3' UTRs using the Interval Tree
#. Remove any reads with alignments to more than one single 3' UTR
#. Remove any additional alignments to non-UTR regions

.. image:: img/multimappers.png
   :width: 800px
   
.. _sample-file:
   
Sample file format
------------------

Sample description files are a convenient way to provide multiple sample files as well as corresponding sample information to *slamdunk*. A sample description
file is a tab or comma separated plain text file with the following columns:

================== ========  ===============================================================================================
Column             Datatype  Description
================== ========  ===============================================================================================
Filepath           String    Path to raw reads in BAM/fasta(gz)/fastq(gz) format
Sample description String    Description of the sample
Sample type        String    Type of sample. Only relevant for timecourses (pulse/chase) otherwise use any value 
Timepoint          Integer   Timepoint of the sample in minutes. Again only relevant for timecourses otherwise use any value
================== ========  ===============================================================================================

.. _tcount-file:

Tcount file format
------------------

The *tcount* file is the central output file of *slamdunk*. It contains all results, conversion rates and other statistics for each UTR which is the
main starting point for any subsequent analysis that will follow e.g. transcript half-life estimates or DE analysis.

*Tcount* files are essentially tab-separated text files containing one line entry per 3' UTR supplied by the user. Each entry contains the following
columns:

===============  ========  ===================================================================================
Column           Datatype  Description
===============  ========  ===================================================================================
chromosome       String    Chromosome on which the 3' UTR resides
start            Integer   Start position of the 3' UTR
end              Integer   End position of the 3' UTR
name             String    Name or ID of the 3' UTR supplied by the user
length           Integer   Length of the 3' UTR
strand           String    Strand of the 3' UTR
conversionRate   Float     Average conversion rate (see :ref:`conversionRates`)
readsCPM         Float     Normalized number of reads as counts per million
Tcontent         Integer   Number of Ts in the 3' UTR (As for - strand UTRs)
coverageOnTs     Integer   Cumulative coverage on all Ts in the 3' UTR
conversionsOnTs  Integer   Cumulative number of T>C conversions in the 3' UTR
readCount        Integer   Number of reads aligning to the 3' UTR
tcReadCount      Integer   Number of reads with T>C conversions aligning to the 3' UTR
multimapCount    Integer   Number of reads considered as multimappers aligning to the 3' UTR
===============  ========  ===================================================================================

Here is an example:

.. code:: bash

   chr13   14197734        14199362        ENSMUSG00000039219      1628    +       0.0466045272969 4.00770947448   587     751     35      59      29      0
   chr2    53217389        53219220        ENSMUSG00000026960      1831    +       0.0270802560315 28.3936027175   709     6093    165     418     138     2
   chr2    53217389        53218446        ENSMUSG00000026960      1057    +       0.0268910814471 23.9783295677   407     5169    139     353     118     0
   chr3    95495394        95495567        ENSMUSG00000015522      173     +       0.0290697674419 1.08683646766   53      172     5       16      5       0
   chr3    95495394        95497237        ENSMUSG00000015522      1843    +       0.0253636702723 11.6834920273   584     2681    68      172     55      4
   chr6    113388777       113389056       ENSMUSG00000079426      279     +       0.0168514412417 16.7780379694   71      2255    38      247     38      3


.. _cB-file:

cB file format
------------------

The *cB* file is a new, optional output of SLAMDUNK, introduced in version 0.5.0. cB stands for "counts Binomial", and is a tidy table that 
is designed to support mixture modeling, a statistically rigorous strategy for estimating the fraction of reads from a given UTR that were
from metabolically labeled reads. This analysis strategy was originally proposed in `Schofield et al., 2018 <https://www.nature.com/articles/nmeth.4582>`_ 
and implemented in software like `GRAND-SLAM <https://academic.oup.com/bioinformatics/article/34/13/i218/5045735?login=true>`_ and 
later `bakR <https://rnajournal.cshlp.org/content/29/7/958.abstract>`_. Mixture modeling overcomes the limitations of using a single T>C conversion
cutoff to classify reads as labeled vs. unlabeled (e.g., RT/sequencing errors in reads from unlabeled RNA, low metabolic label incorporation rates,
etc.). bakR can be provided a *cB* file as input to perform mixture modeling for you. 

*cB* files are essentially comma-separated text files containing one line entry per group of reads with identical "information content".
Information content refers to the UTR from which the read originated, as well as the mutational (T>C) and nucleotide content (T) of the read. Thus,
the columns contained in this file are as follows:

===============  ========  ===================================================================================
Column           Datatype  Description
===============  ========  ===================================================================================
chromosome       String    Chromosome on which the 3' UTR resides
start            Integer   Start position of the 3' UTR
end              Integer   End position of the 3' UTR
name             String    Name or ID of the 3' UTR supplied by the user
strand           String    Strand of the 3' UTR
TC               Integer   Number of T>C conversions in the read
nT               Integer   Number of reference Ts covered by the read
n                Integer   Number of reads that share all of the information described in the other columns
===============  ========  ===================================================================================

Here is an example:

.. code:: bash

   chr13   14197734        14199362        ENSMUSG00000039219      +       0  25 10
   chr13   14197734        14199362        ENSMUSG00000039219      +       0  26 28
   chr13   14197734        14199362        ENSMUSG00000039219      +       0  28 5
   chr13   14197734        14199362        ENSMUSG00000039219      +       1  20 3
   chr13   14197734        14199362        ENSMUSG00000039219      +       1  25 15
   chr6    113388777       113389056       ENSMUSG00000079426      +       0  30 2