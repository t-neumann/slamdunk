
bin/slamdunk all -r slamdunk/test/data/ref.fa -b slamdunk/test/data/actb.bed -o slamdunk/test/data/output -rl 100 -mbq 27 -5 0 slamdunk/test/data/reads.fq

diff slamdunk/test/data/reads_slamdunk_mapped_filtered_tcount.tsv slamdunk/test/data/output/count/reads_slamdunk_mapped_filtered_tcount.tsv 

