Changelog
=========

**Version 0.3.0**

*Changes:*

* Travis-CI testing added.
* Automated versionable Docker container builds.
* Static 2 T>C conversion cutoff multiTCSTringency parameter replaced by dynamic conversion-threshold parameter
* Positional track module added to Alleyoop to produce genome-wide positional T>C conversion rate bigWigtracks
* Read-separator module added to Alleyoop to separated T>C reads from non T>C reads into bam-files.
* Docker base image set to ubuntu:18.04

*Bugfixes:*

* Supplementary alignment flags unset in filtered bam
* Fixed numpy version removed to satisfy pandas dependency
* Minimum baseQ filter now inclusive which is the expected behaviour
* Reference check for short chromosomes introduced
* Backwards compatibility to pysam 0.8.3 introduced

**Version 0.2.4**

*Changes:*

* License updated to AGPLv3
* Deduplicator now only marks duplicates instead of removing them. Useful for deduplication QA like `dupRadar <https://bioconductor.org/packages/release/bioc/html/dupRadar.html>`_. Allows for deduplication of T>C-fraction of reads only.
* `Singularity <http://singularity.lbl.gov/>`_ build-file branch added.

*Bugfixes:*

* readInRegion fixed: Switched from fetch(region="chr:start-end" (1-based) to fetch(reference=chr, start=start, end=end) (0-based). Intervals close to the ends of a chromosome are now correctly handeled and not skipped.
* Hotfix to fix executable permissions on some R plotting scripts
* Fixed alleyoop collapse to ignore all MLE columns for now
* R package dependencies fixed for IMP cluster
* include_package_data removed for proper contrib import

**Version 0.2.3**

*Bugfixes:*

* Hotfix to include README.md in sdist packaging
* Hotfix to fix executable permissions on some R plotting scripts
* Fixed alleyoop collapse to ignore all MLE columns for now

**Version 0.2.2**

*Changes:*

* Slamdunk supports now `MultiQC <http://multiqc.info/>`_ report creation
* Default set to local mapping
* Base-quality cutoff implemented throughout the entire *slamdunk* / *alleyoop* workflow
* BAQ computation disabled for samtools mpileup
* Dedup now supports selection of reads containing a defined # T>C conversions
* All documentation moved to gh-pages branch for slicker repo
* CPMs now normalized to filtered reads
* Docker image now build automatically to head revision
* Official webpage introduced instead of readthedocs.org hosting
* Evaluation module for count files added to splash

*Bugfixes:*

* PCAs for single sample working now
* Dedup working again
* Sample files may contain empty lines now


**Version 0.2.1**

*Changes:*

* PCA will be added to alleyoop summary if count folder is provided 
* Reads in UTRs will be added to alleyoop summary if count folder is provided
* Comment tags added to summary, rates, utrrates, tcperreadpos, tcperutrpos to be recognized by `MultiQC <http://multiqc.info/>`
* NextGenMap udated to version 0.5.2
* Slamdunk now checks NextGenMap version
* Simulation module (splash) added: simulates slamdunk samples from a set of UTRs (BED file) and evalutes results computed by SlamDunk. See `splash --help` for more information 
* `-i/--sample-index` is now 1-based instead of 0-based
* Columns for conversion rate confidence interval added to count files (not used yet)
* Filter summary integrated into log
* Added sample names now in globalRatePlotter

*Bugfixes:*

* Number of sequenced reads fixed for BAM files computed with `-ss/--skip-sam` option
* NextGenMap memory leak fixed

**Version 0.2.0** 

Due to major changes v0.2.0 is not compatible with v0.1.0. 
Please rerun alle your samples with slamdunk 0.2.0!

*Changes:*

* Option to supply sample info via samplesheet instead of plain raw reads added. 
* Added @RG tag to bam files containing sample info for best-practice and easier summary calculation. 
* Flagstat files removed. 
* --skip-sam flag added to only output bam. 
* --sample-index flag implemented for cluster distribution. 
* Multi-TC stringency now also available for utr rate plots. 
* \*stats prefixes removed from alleyoop modules.
* Alleyoop summary module implemented. Summary contains sample infos (if specified via samplesheet) and number of sequenced, mapped and filtered reads.  
* Alleyoop merge module allows now merging upon column names. 
* Sample info and reference + checksum now documented in tcount files. 
* Version.py added. 

*Bugfixes*:
 
* NGM interspersed order of multimappers fixed. 
* Random seed from filter module deleted. 
* Alleyoop rates display error of bars higher than ylim fixed. 
* Min base quality now again propagated to all modules. 
* Auto-scaling error fixed in globalRatePlotter.  

**Version 0.1.0** 

* Initial pre-release.
