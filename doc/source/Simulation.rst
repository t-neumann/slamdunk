Simulating SlamSeq data
=====

Example
^^^^^^^

.. code:: bash

    #!bash
    set -x
	F=simulation_2/
	B=finalAnnotation_test_cut_chrM_correct_100
	RL=50
	COV=100
	SN=pooja_UTR_annotation_examples_sample
	
	#rm -rf $F
	#~/bin/slamdunk/bin/slamsim preparebed -b ${B}.bed -o ${F} -l ${RL}
	#~/bin/slamdunk/bin/slamsim utrs -b ${F}/${B}_original.bed -o $F -r /project/ngs/philipp/slamseq/ref/GRCm38.fa
	
	#rm -rf ${F}/slamdunk ${F}/*sample*
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 0 --sample-name ${SN}_1_0min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 15 --sample-name ${SN}_2_15min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 30 --sample-name ${SN}_3_30min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 60 --sample-name ${SN}_4_60min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 180 --sample-name ${SN}_5_180min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 360 --sample-name ${SN}_6_360min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 720 --sample-name ${SN}_7_720min
	#~/bin/slamdunk/bin/slamsim reads -b ${F}/${B}_original_utrs.bed -l ${RL} -o $F -cov ${COV} -t 1440 --sample-name ${SN}_8_1440min
	
	#~/bin/slamdunk/bin/slamdunk all -r /project/ngs/philipp/slamseq/ref/GRCm38.fa -b ${F}/${B}_original_utrs.bed -rl 55 -o ${F}/slamdunk ${F}/*.bam
	
	~/bin/slamdunk/bin/slamsim plot.conversions -sim ${F} -slam ${F}/slamdunk/count/ -o ${F}/eval/conversion_rate_eval_plots.pdf
	
	~/bin/slamdunk/bin/slamsim plot.halflifes -sim ${F} -slam ${F}/slamdunk/count/ -t 0,15,30,60,180,360,720,1440 -o ${F}/eval/halflife_per_gene_eval_plots.pdf -b ${F}/${B}_original_utrs.bed
	

