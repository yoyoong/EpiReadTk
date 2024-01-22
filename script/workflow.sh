echo "biscuit epiread -B 5.snp_bed/snp.bed -N 0.refer_fastq/hg38.fa 3.bam/Z000000PX_sorted.bam | sort -k1,1 -k2,2n > 6.epibed/Z000000PX.epibed" \
| qsub -q b1.q@fnode006.sibcb.ac.cn -N biscuit_Z000000PX -V -cwd -o /sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/log/biscuit_Z000000PX.log -j y

