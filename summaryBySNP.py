import os
from pathlib import Path
from common.Namespace import Namespace
from common.Region import Region
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile
from common.FastaFile import FastaFile
from common.BedFile import BedFile
from common.Constant import CODE, CPG_DICT, SNP_DICT, SNP_DICT_INV, SNP_REPLACE_ATCG_DICT
from common.OutputFile import OutputFile
from common.Util import calculate_MHL, calculate_MBS, calculate_Entropy, calculate_R2
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, find
import datetime
import gzip
from tqdm.autonotebook import tqdm

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
        self.k4Plus = args.k4Plus
        self.cpgCov = args.cpgCov
        self.r2Cov = args.r2Cov

    def check_args(self):
        assert isinstance(self.maxK, int), 'maxK should be int'
        assert isinstance(self.K, int), 'K should be int'
        assert isinstance(self.minK, int), 'minK should be int'
        assert Path(os.path.dirname(self.outputDir)).exists(), 'output directory does not exist'
        assert 1 <= self.minK <= 10, 'minK should be in 1 to 10'
        assert 1 <= self.maxK <= 10, 'maxK should be in 1 to 10'
        assert self.maxK > self.minK, 'maxK should be larger than minK'
        # assert 3 <= self.K <= 5, 'K：the default is 4 and values must be between 3 and 5'
        assert self.region or self.bedFile, 'you should input bedPath or region'
        assert (not self.region) or (not self.bedFile), 'you should only input bedPath or region'
        assert self.strand == 'both' or self.strand == 'plus' or self.strand == 'minus', 'strand should be both, plus or minus'

    def summary_by_region(self, region):
        self.epibed_info = self.epibedFile.query_by_region(region)
        self.cpg_snp_position = self.epibedFile.get_cpg_snp_position()
        self.cpg_snp_matrix, self.strand_list = self.epibedFile.build_sparse_cpg_snp_matrix()
        cpg_pos_list = self.cpgFile.query_by_region(self.region)

        # get the position has and only has 2 kind ATCG replacement snp
        cpg_snp_matrix_row, cpg_snp_matrix_col, cpg_snp_matrix_data = find(self.cpg_snp_matrix)
        atcg_mask = (cpg_snp_matrix_data <= CODE.G_REPLACE.value) & (cpg_snp_matrix_data >= CODE.A_REPLACE.value)
        atcg_replace_matrix = coo_matrix(
            (cpg_snp_matrix_data[atcg_mask], (cpg_snp_matrix_row[atcg_mask], cpg_snp_matrix_col[atcg_mask])),
            shape=self.cpg_snp_matrix.shape)
        nonzero_row, nonzero_col = atcg_replace_matrix.nonzero()
        atcg_replace_index = np.unique(nonzero_col)
        # insert和delete暂不处理

        # cpg information dictionatry grouped by ATCG replacement and indexed position
        pos_snp_cpg_dict = dict()
        cpg_mask = (cpg_snp_matrix_data == CODE.UNMETHYLATED.value) | (cpg_snp_matrix_data == CODE.METHYLATED.value)
        cpg_info_matrix = coo_matrix(
            (cpg_snp_matrix_data[cpg_mask], (cpg_snp_matrix_row[cpg_mask], cpg_snp_matrix_col[cpg_mask])),
            shape=self.cpg_snp_matrix.shape)
        for i in tqdm(atcg_replace_index, desc='Get cpg info'):
            atcg_replace_info = atcg_replace_matrix.getcol(i)
            if len(np.unique(atcg_replace_info.data)) == 2:
                snp_cpg_dict = dict()
                for snp_code in atcg_replace_info.data:
                    # if position in cpg_pos_list: # 如果即是cpg位点也是snp位点，该位点的cpg信息忽略
                    cpg_info_index = self.cpg_snp_matrix.getcol(i).multiply(self.cpg_snp_matrix.getcol(i) == snp_code).nonzero()[0]
                    cpg_info = cpg_info_matrix.todense()[cpg_info_index]
                    if not np.all(cpg_info == 0):
                        snp_cpg_dict[snp_code] = cpg_info

                position = self.cpg_snp_position[i]
                if len(snp_cpg_dict) == 2:
                    pos_snp_cpg_dict[position] = snp_cpg_dict
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " 4!")

        output_file = OutputFile(os.path.join(self.outputDir, self.tag + ".txt"))
        head_base_list = ["chrom", "pos", "refer", "real", "nReads", "mBase", "cBase", "tBase", "K4plus", "nDR", "nMR"]
        metrics_list = self.metrics.split(" ")
        column_list = head_base_list + metrics_list
        output_file.write_head(column_list)
        for pos, snp_cpg_dict in tqdm(sorted(pos_snp_cpg_dict.items()), desc='Calculate'):
            for snp_code, cpg_info in sorted(snp_cpg_dict.items()):
                line_str = self.calculate_by_site(pos, snp_code, cpg_info, metrics_list)
                output_file.write_line(line_str)

        # print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " 5!")

    def calculate_by_site(self, site_pos:str, snp_code:int, cpg_info: np.ndarray, metrics_list:list):
        chrom = self.region.chr
        pos = site_pos
        refer = self.fastaFile.query_by_region(self.region.chr, site_pos, site_pos + 1).upper()
        real = SNP_DICT_INV[snp_code]
        nReads = 0 # 总read个数
        mBase = 0 # 甲基化位点个数
        cBase = 0 # 存在甲基化的read中的未甲基化位点个数
        tBase = 0 # 总位点个数
        K4plus = 0 # 长度大于等于K个位点的read个数
        nDR = 0 # 长度大于等于K个位点且同时含有甲基化和未甲基化位点的read个数
        nMR = 0 # 长度大于等于K个位点且含有甲基化位点的read个数
        for read in cpg_info:
            nReads += 1
            mBase += np.count_nonzero(read == CODE.METHYLATED.value)
            tBase += np.count_nonzero(read != 0)
            if CODE.METHYLATED.value in read:
                cBase += np.count_nonzero(read == CODE.UNMETHYLATED.value)
            if np.count_nonzero(read != 0) > self.K:
                K4plus += 1
                if CODE.METHYLATED.value in read:
                    nMR += 1
                    if CODE.UNMETHYLATED.value in read:
                        nDR += 1

        line_str = f'{chrom}\t{pos}\t{refer}\t{real}\t{nReads}\t{mBase}\t{cBase}\t{tBase}\t{K4plus}\t{nDR}\t{nMR}\t'
        for metrics in metrics_list:
            if metrics == "MM":
                if nReads < self.cpgCov:
                    line_str += str(float('nan')) + '\t'
                else:
                    MM = "{:.4f}".format(mBase / tBase)
                    line_str += f'{MM}\t'
            elif metrics == "PDR":
                if K4plus < self.k4Plus:
                    line_str += str(float('nan')) + '\t'
                else:
                    PDR = "{:.4f}".format(nDR / K4plus)
                    line_str += f'{PDR}\t'
            elif metrics == "CHALM":
                if K4plus < self.k4Plus:
                    line_str += str(float('nan')) + '\t'
                else:
                    CHALM = "{:.4f}".format(nMR / K4plus)
                    line_str += f'{CHALM}\t'
            elif metrics == "MHL":
                if K4plus < self.k4Plus:
                    line_str += str(float('nan')) + '\t'
                else:
                    MHL = "{:.4f}".format(calculate_MHL(self.minK, self.maxK, cpg_info))
                    line_str += f'{MHL}\t'
            elif metrics == "MCR":
                if nReads < self.cpgCov:
                    line_str += str(float('nan')) + '\t'
                else:
                    MCR = "{:.4f}".format(cBase / tBase)
                    line_str += f'{MCR}\t'
            elif metrics == "MBS":
                if K4plus < self.k4Plus:
                    line_str += str(float('nan')) + '\t'
                else:
                    MBS = "{:.4f}".format(calculate_MBS(self.K, cpg_info))
                    line_str += f'{MBS}\t'
            elif metrics == "Entropy":
                if K4plus < self.k4Plus:
                    line_str += str(float('nan')) + '\t'
                else:
                    Entropy = "{:.4f}".format(calculate_Entropy(self.K, cpg_info))
                    line_str += f'{Entropy}\t'
            elif metrics == "R2":
                R2 = "{:.4f}".format(calculate_R2(cpg_info))
                line_str += f'{R2}\t'

        return line_str.rstrip() + '\n'

    def calculate_all(self, ):
        cpg_matrix = pd.read_csv(self.cpg_path, sep='\t', header=None).iloc[:, :2]


def main(args):
    print("Run summaryBySNP start!")
    summaryBySNP = SummaryBySNP(args)
    summaryBySNP.check_args()

    # get region list from input args
    if summaryBySNP.region is not None:
        summaryBySNP.summary_by_region(summaryBySNP.region)
    elif summaryBySNP.bedFile is not None:
        region_list = summaryBySNP.bedFile.get_region_list()
        for region in region_list:
            summaryBySNP.summary_by_region(region)
    else:
        summaryBySNP.calculate_all()

    print("Run summaryBySNP end!")

if __name__ == '__main__':
    args = Namespace()
    args.epibedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args.epibedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args.cpgPath = "/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg19/hg19_CpG.gz"
    args.fastaPath = "/sibcb2/bioinformatics/iGenome/Bismark/hg19/hg19.fa"
    args.region = "chr1:923100-10924200"
    # args.bedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/bed_file/test.bed"
    args.outputDir = "/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir"
    args.tag = "summaryBySNP.test"
    args.metrics = "MM PDR CHALM MHL MCR MBS Entropy"
    # args.metrics = "Entropy"
    args.minK = 1
    args.maxK = 10
    args.K = 4
    args.cutReads = True
    args.strand = "both"
    args.k4Plus = 4
    args.cpgCov = 1
    args.r2Cov = 1
    main(args)
