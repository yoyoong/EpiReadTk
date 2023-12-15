import argparse
import os
import warnings
warnings.filterwarnings('ignore')

import gzip
from main import HelpMessage
from common.OutputFile import OutputFile
from tqdm.autonotebook import tqdm

# Argument parser
parser = argparse.ArgumentParser(description='Get snp bed file from vcf file')
parser.add_argument('--vcf_path', type=str, required=False, help=HelpMessage.vcf_path,
                    default='/sibcb2/bioinformatics/iGenome/GATK/bundle/hg38/dbsnp_146.hg38.vcf.gz')
parser.add_argument('--snp_id', type=str, required=False, help=HelpMessage.snp_id,
                    default='/sibcb2/bioinformatics/snpWorkflow/out/common_snp_id.txt')
parser.add_argument('--outputDir', type=str, required=False, help=HelpMessage.outputDir,
                    default='/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/5.snp_bed')
parser.add_argument('--tag', type=str, required=False, help=HelpMessage.tag + " [snp]", default='snp')
args = parser.parse_args()

vcf_path = args.vcf_path
snp_id_path = args.snp_id
outputDir = args.outputDir
tag = args.tag

# get snp id
with open(snp_id_path, 'r') as file:
    snp_id_lines = file.readlines()
snp_id_list = set([line.strip() for line in snp_id_lines])

# get snp reference file and generate snp bed file
snp_bed_path = os.path.join(outputDir, f"{tag}.bed")
snpFile = OutputFile(snp_bed_path)
# column_list = ["chrom", "start", "end", "refer_base", "alter_base", "genotype", "snp_id", "coverage", "frequency"]
# snpFile.write_head(column_list)
with tqdm(gzip.open(vcf_path, 'rt'), desc="Read vcf file") as file:
    for line in file:
        if line.startswith("#"):
            continue
        else:
            info = line.split("\t")
            snp_id = info[2]
            if snp_id not in snp_id_list:
                continue
            chrom = info[0]
            start = int(info[1]) - 1 # 0-based start (BED format)
            end = start + 1 # 1-based non-inclusive end (BED format)
            refer_base = info[3]
            alter_base = info[4]
            genotype = "0/1"
            coverage = "."
            frequency = "0.5"
            snp_line = f'{chrom}\t{start}\t{end}\t{refer_base}\t{alter_base}\t{genotype}\t{snp_id}\t{coverage}\t{frequency}\n'
            snpFile.write_line(snp_line)
