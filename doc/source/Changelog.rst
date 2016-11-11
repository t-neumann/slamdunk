Changelog
=========

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
* Min base quality now again propagated to all modules. Auto-scaling error fixed in globalRatePlotter.  

**0.1.0:** 

* Initial pre-release.