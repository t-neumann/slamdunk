Changelog
=========

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