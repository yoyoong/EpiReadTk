import os
from pathlib import Path
from common.Namespace import Namespace
from common.Region import Region
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile
from common.FastaFile import FastaFile
from common.BedFile import BedFile
from common.Constant import CODE, CPG_DICT, SNP_DICT, SNP_DICT_INV, SNP_REPLACE_ATCG_DICT
import numpy as np
import pandas as pd
import copy

SHIFT = 500

class SummaryBySNP:
    def __init__(self, args):
        self.epibedFile = EpibedFile(args.epibedPath)
        self.cpgFile = CpGFile(args.cpgPath, SHIFT)
        self.fastaFile = FastaFile(args.fastaPath)
        self.region = Region(args.region) if args.__contains__("region") else None
        self.bedFile = BedFile(args.bedPath, 0) if args.__contains__("bedPath") else None
        self.outputDir = args.outputDir
        self.tag = args.tag
        self.metrics = args.metrics
        self.minK = args.minK
        self.maxK = args.maxK
        self.K = args.K
        self.cutReads = args.cutReads
        self.strand = args.strand
        self.r2Cov = args.r2Cov

    def check_args(self):
        assert isinstance(self.maxK, int), 'maxK should be int'
        assert isinstance(self.K, int), 'K should be int'
        assert isinstance(self.minK, int), 'minK should be int'
        assert Path(os.path.dirname(self.outputDir)).exists(), 'output directory does not exist'
        assert 1 <= self.minK <= 10, 'minK should be in 1 to 10'
        assert 1 <= self.maxK <= 10, 'maxK should be in 1 to 10'
        assert self.maxK > self.minK, 'maxK should be larger than minK'
        assert 3 <= self.K <= 5, 'K：the default is 4 and values must be between 3 and 5'
        assert self.region or self.bedFile, 'you should input bedPath or region'
        assert (not self.region) or (not self.bedFile), 'you should only input bedPath or region'
        assert self.strand == 'both' or self.strand == 'plus' or self.strand == 'minus', 'strand should be both, plus or minus'

    def summary_region(self, region):
        self.epibed_info = self.epibedFile.query_by_region(region)
        self.cpg_snp_position = self.epibedFile.get_cpg_snp_position()
        self.cpg_snp_matrix, self.strand_list = self.epibedFile.build_cpg_snp_matrix()
        cpg_pos_list = self.cpgFile.query_by_region(self.region)
        pos_num = self.cpg_snp_matrix.shape[0]
        read_num = self.cpg_snp_matrix.shape[1]

        # get the position has and only has 2 kind ATCG replacement snp
        atcg_replace_mask = np.isin(self.cpg_snp_matrix, list(SNP_REPLACE_ATCG_DICT.values()))
        atcg_replace_matrix = self.cpg_snp_matrix * atcg_replace_mask
        position_atcg_replace_info = [np.unique(col[col != 0]) for col in atcg_replace_matrix.T]
        position_atcg_replace_num = [len(item) for item in position_atcg_replace_info]
        position_for_index = np.array(self.cpg_snp_position)[np.array(position_atcg_replace_num) == 2]
        # insert和delete是否要处理？

        # cpg information dictionatry grouped by ATCG replacement and indexed position
        pos_snp_cpg_dict = dict()
        for i in range(len(position_atcg_replace_info)):
            if len(position_atcg_replace_info[i]) == 2:
                snp_cpg_dict = dict()
                for snp_code in position_atcg_replace_info[i]:
                    # if position in cpg_pos_list: # 如果即是cpg位点也是snp位点，该位点的cpg信息是忽略还是作为未甲基化处理？
                    cpg_mask = np.isin(self.cpg_snp_matrix, list(CPG_DICT.values()))
                    cpg_info = (self.cpg_snp_matrix * cpg_mask)[self.cpg_snp_matrix.T[i] == snp_code, :]
                    if not np.all(cpg_info == 0):
                        snp_cpg_dict[snp_code] = cpg_info

                position = self.cpg_snp_position[i]
                if len(snp_cpg_dict) > 0:
                    pos_snp_cpg_dict[position] = snp_cpg_dict

        print("pos_snp_cpg_dict")

def main(args):
    print("Run summaryBySNP start!")
    summaryBySNP = SummaryBySNP(args)
    summaryBySNP.check_args()

    # get region list from input args
    region_list = list()
    if summaryBySNP.bedFile is not None:
        region_list = summaryBySNP.bedFile.get_region_list()
    else:
        region_list.append(summaryBySNP.region)

    for region in region_list:
        summaryBySNP.summary_region(region)

    print("Run summaryBySNP end!")

if __name__ == '__main__':
    args = Namespace()
    args.epibedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args.epibedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args.cpgPath = "/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg19/hg19_CpG.gz"
    args.fastaPath = "/sibcb2/bioinformatics/iGenome/Bismark/hg19/hg19.fa"
    args.region = "chr1:230177036-230178500"
    # args.bedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/bed_file/test.bed"
    args.outputDir = "/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/summaryBySNP_outputDir"
    args.tag = "summaryBySNP.test"
    args.metrics = "MM"
    args.minK = 1
    args.maxK = 10
    args.K = 4
    args.cutReads = True
    args.strand = "both"
    args.r2Cov = 20
    main(args)
