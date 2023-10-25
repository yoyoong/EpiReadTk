import os
from pathlib import Path
from common.Namespace import Namespace
from common.Region import Region
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile
from common.FastaFile import FastaFile
from common.BedFile import BedFile
from common.Constant import CODE, CPG_DICT, SNP_DICT, SNP_DICT_INV, SNP_INSERT_DICT, SNP_REPLACE_ATCG_DICT
from common.OutputFile import OutputFile
from common.Util import calculate_MHL, calculate_MBS, calculate_Entropy, calculate_R2
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, find
import datetime
import gzip
from tqdm.autonotebook import tqdm
import itertools
import sys

class SummaryBySNP:
    def __init__(self, args):
        self.epibedFile = EpibedFile(args.epibedPath)
        self.cpgFile = CpGFile(args.cpgPath)
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
        # assert self.region or self.bedFile, 'you should input bedPath or region'
        assert (not self.region) or (not self.bedFile), 'you should only input bedPath or region'
        assert self.strand == 'both' or self.strand == 'plus' or self.strand == 'minus', 'strand should be both, plus or minus'

    def summary_by_region(self, region):
        self.epibed_info = self.epibedFile.query_by_region(region)
        self.cpg_snp_position = self.epibedFile.get_cpg_snp_position()
        self.cpg_snp_matrix = self.epibedFile.build_sparse_cpg_snp_matrix()
        # cpg_pos_list = self.cpgFile.query_by_region(self.region)
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + " 2!")

        # get the position has ATCG replacement snp
        cpg_snp_matrix_row, cpg_snp_matrix_col, cpg_snp_matrix_data = find(self.cpg_snp_matrix)
        atcg_mask = (cpg_snp_matrix_data <= CODE.G_REPLACE.value) & (cpg_snp_matrix_data >= CODE.A_REPLACE.value)
        atcg_replace_matrix = coo_matrix(
            (cpg_snp_matrix_data[atcg_mask], (cpg_snp_matrix_row[atcg_mask], cpg_snp_matrix_col[atcg_mask])),
            shape=self.cpg_snp_matrix.shape)
        atcg_replace_cols = np.unique(atcg_replace_matrix.nonzero()[1]) # column index include ATCG replacement
        # insert和delete暂不处理
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + " 3!")

        # cpg information dictionatry grouped by ATCG replacement and indexed position
        pos_snp_cpg_dict = dict()
        cpg_mask = (cpg_snp_matrix_data == CODE.UNMETHYLATED.value) | (cpg_snp_matrix_data == CODE.METHYLATED.value)
        cpg_info_matrix = coo_matrix(
            (cpg_snp_matrix_data[cpg_mask], (cpg_snp_matrix_row[cpg_mask], cpg_snp_matrix_col[cpg_mask])),
            shape=self.cpg_snp_matrix.shape)
        cpg_info_matrix = cpg_info_matrix.todense()
        for i in tqdm(atcg_replace_cols, desc='Get cpg info'):
            atcg_replace_info = atcg_replace_matrix.getcol(i)
            if len(np.unique(atcg_replace_info.data)) == 2:
                snp_cpg_dict = dict()
                cpg_snp_matrix_col = self.cpg_snp_matrix.getcol(i)
                for snp_code in atcg_replace_info.data:
                    # if position in cpg_pos_list: # 如果即是cpg位点也是snp位点，该位点的cpg信息忽略
                    cpg_info_index = cpg_snp_matrix_col.multiply(cpg_snp_matrix_col == snp_code).nonzero()[0]
                    cpg_info = cpg_info_matrix[cpg_info_index]
                    if not np.all(cpg_info == 0):
                        snp_cpg_dict[snp_code] = np.array(cpg_info)

                position = self.cpg_snp_position[i]
                if len(snp_cpg_dict) == 2:
                    pos_snp_cpg_dict[position] = snp_cpg_dict
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + " 4!")

        output_file = OutputFile(os.path.join(self.outputDir, self.tag + ".txt"))
        head_base_list = ["chrom", "pos", "refer", "real", "nReads", "mBase", "cBase", "tBase", "K4plus", "nDR", "nMR"]
        metrics_list = self.metrics.split(" ")
        column_list = head_base_list + metrics_list
        output_file.write_head(column_list)
        for pos, snp_cpg_dict in tqdm(sorted(pos_snp_cpg_dict.items()), desc='Calculate'):
            for snp_code, cpg_info in sorted(snp_cpg_dict.items()):
                line_str = self.calculate_by_site(pos, snp_code, cpg_info, metrics_list)
                output_file.write_line(line_str)

        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + " 5!")

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
            if np.count_nonzero(read != 0) >= self.K:
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

    def separate_rle_string(self, rle_string: str):
        rle_list = []
        for key, group in itertools.groupby(rle_string, str.isdigit):
            value = "".join(group)
            if value.isdigit():
                rle_list.append(value)
            else:
                for char in value:
                    rle_list.append(char)
        return rle_list

    def calculate_all(self):
        head_base_list = ["chrom", "pos", "refer", "real", "nReads", "mBase", "cBase", "tBase", "K4plus", "nDR", "nMR"]
        metrics_list = self.metrics.split(" ")
        if "MHL" in metrics_list:
            for i in range(self.minK - 1, self.maxK):
                head_base_list.append(f'methKmers{i}')
                head_base_list.append(f'totalKmers{i}')
        if "MBS" in metrics_list:
            head_base_list.append('mbsNum')
        if "Entropy" in metrics_list:
            for i in range(args.K):
                for j in range(args.K):
                    head_base_list.append(f'kmer{i}_{j}')
            head_base_list.append('kmerAll')

        column_list = head_base_list + metrics_list
        raw_data = pd.DataFrame(columns=column_list)

        for read in tqdm(gzip.open(self.epibedFile.epibed_path, 'rb'), "Read EpiRead"):
            line = read.decode().split('\t')
            start_pos = int(line[1]) + 1  # 0-base to 1-base
            read_strand = line[5]
            if self.strand != 'both':
                if self.strand != read_strand:
                    continue

            # get the snp number in this line
            snp_num = 0
            for s in line[8]:
                if s in SNP_REPLACE_ATCG_DICT.keys(): # 先考虑ATCG replacement
                    snp_num += 1
            if snp_num == 0:
                continue

            if int(line[1]) < 882400:
                continue
            if int(line[1]) > 1000000:
                break;

            cpg_info = self.separate_rle_string(line[6])
            snp_info = self.separate_rle_string(line[8])

            # get the snp position and label list
            snp_pos_list = []
            snp_label_list = []
            snp_move_length = 0
            snp_insert_length = 0
            for i in range(len(snp_info)):
                if snp_info[i].isalpha():
                    if snp_info[i] in SNP_DICT.keys():
                        if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                            for s in range(int(snp_info[i + 1])):
                                position = start_pos + snp_move_length - snp_insert_length + s
                                snp_pos_list.append(position)
                                snp_label_list.append(snp_info[i])
                            snp_move_length += int(snp_info[i + 1])
                            if snp_info[i] in SNP_INSERT_DICT.keys():
                                snp_insert_length += int(snp_info[i + 1])
                        else:
                            position = start_pos + snp_move_length - snp_insert_length
                            snp_pos_list.append(position)
                            snp_label_list.append(snp_info[i])
                            snp_move_length += 1
                            if snp_info[i] in SNP_INSERT_DICT.keys():
                                snp_insert_length += 1
                    else:
                        if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                            snp_move_length += int(snp_info[i + 1])
                        else:
                            snp_move_length += 1

            # get the cpg position and label list
            cpg_pos_list = []
            cpg_raw_label_list = []
            cpg_move_length = 0
            cpg_insert_length = 0
            for i in range(len(cpg_info)):
                if cpg_info[i].isalpha():
                    if cpg_info[i] in CPG_DICT.keys():
                        if read_strand == "+":  # reads from OT/CTOT(+) strands, methylation site is in C→T substitution
                            position = start_pos + cpg_move_length - cpg_insert_length
                        else:  # reads from OB/CTOB (-) strands, methylation site is in G→A substitution
                            position = start_pos + cpg_move_length - cpg_insert_length - 1
                        if position not in snp_pos_list:  # 在cpg位点上的snp不考虑
                            cpg_pos_list.append(position)
                            cpg_raw_label_list.append(cpg_info[i])

                    if i + 1 < len(cpg_info) and cpg_info[i + 1].isdigit():
                        cpg_move_length += int(cpg_info[i + 1])
                        if cpg_info[i] in SNP_INSERT_DICT.keys():
                            cpg_insert_length += int(cpg_info[i + 1])
                    else:
                        cpg_move_length += 1
                        if cpg_info[i] in SNP_INSERT_DICT.keys():
                            cpg_insert_length += 1
            if len(cpg_raw_label_list) < 1:
                continue

            # get the cpg label string, '1'-methylated, '0'-unmethylated, '-' is not cover
            region = Region.init(line[0], cpg_pos_list[0], cpg_pos_list[-1] + 1)
            cpg_full_pos_list = self.cpgFile.query_by_region(region)
            cpg_label_list = ["-"] * len(cpg_full_pos_list)
            for i in range(len(cpg_pos_list)):
                index = cpg_full_pos_list.index(cpg_pos_list[i])
                if cpg_raw_label_list[i] == CODE.METHYLATED.label:
                    cpg_label_list[index] = "1"
                elif cpg_raw_label_list[i] == CODE.UNMETHYLATED.label:
                    cpg_label_list[index] = "0"
            cpg_label_str = ''.join(cpg_label_list)

            for i in range(len(snp_pos_list)):
                snp_pos = snp_pos_list[i]
                chrom = line[0]
                pos = snp_pos
                real = snp_label_list[i]
                refer = self.fastaFile.query_by_region(chrom, snp_pos, snp_pos + 1).upper()
                nReads = 1  # 总read个数
                mBase = cpg_label_list.count('1')  # 甲基化位点个数
                tBase = len(cpg_raw_label_list) # 总位点个数
                cBase = 0  # 存在甲基化的read中的未甲基化位点个数
                K4plus = 0  # 长度大于等于K个位点的read个数
                nDR = 0  # 长度大于等于K个位点且同时含有甲基化和未甲基化位点的read个数
                nMR = 0  # 长度大于等于K个位点且含有甲基化位点的read个数
                if mBase > 0:
                    cBase = cpg_label_list.count('0')
                if tBase >= self.K:
                    K4plus = 1
                    if mBase > 0:
                        nMR = 1
                        if cBase > 0:
                            nDR = 1

                methKmers_list = [0] * (self.maxK - self.minK + 1)
                totalKmers_list = [0] * (self.maxK - self.minK + 1)
                if "MHL" in metrics_list:
                    if self.minK > len(cpg_label_list):
                        print("Warning: minK is too large for calculate MHL.")
                    maxK = min(self.maxK, len(cpg_label_list))
                    if (self.minK > maxK):
                        print("Error: maxK is too large.\n")
                        sys.exit(-1)
                    for i in range(self.minK - 1, maxK):
                        for j in range(len(cpg_label_list) - i):
                            cpg_label_concat = cpg_label_str[j:(j + i + 1)]
                            if '-' not in cpg_label_concat:
                                cpg_full_meth = '1' * (i + 1)
                                if cpg_label_concat == cpg_full_meth:
                                    methKmers_list[i] += 1
                                totalKmers_list[i] += 1

                mbsNum = 0.0
                if "MBS" in metrics_list:
                    if len(cpg_label_list) >= self.K:
                        cpg_label_split_null = cpg_label_str.split('-')
                        for item in cpg_label_split_null:
                            cpg_label_split = list(filter(None, item.split('0')))
                            mbsNum = sum([len(x) ** 2 for x in cpg_label_split]) / len(cpg_label_list) ** 2

                kmer_list = [0] * (self.K * self.K) # every kmer pattern number
                kmerAll = 0
                if "Entropy" in metrics_list:
                    if len(cpg_label_list) >= self.K:
                        for i in range(len(cpg_label_str) - self.K + 1):
                            kmer = cpg_label_str[i: i + self.K]
                            if '-' not in kmer:
                                kmer_list[int(kmer, 2)] += 1
                                kmerAll += 1

                data_index = chrom + str(pos) + real + refer
                try: # try success mean data_index exist
                    raw_data.loc[data_index, 'nReads'] += nReads
                    raw_data.loc[data_index, 'mBase'] += mBase
                    raw_data.loc[data_index, 'tBase'] += tBase
                    raw_data.loc[data_index, 'cBase'] += cBase
                    raw_data.loc[data_index, 'K4plus'] += K4plus
                    raw_data.loc[data_index, 'nDR'] += nDR
                    raw_data.loc[data_index, 'nMR'] += nMR
                    if "MHL" in metrics_list:
                        for i in range(self.minK - 1, self.maxK):
                            raw_data.loc[data_index, f'methKmers{i}'] += methKmers_list[i]
                            raw_data.loc[data_index, f'totalKmers{i}'] += totalKmers_list[i]

                    if "MBS" in metrics_list:
                        raw_data.loc[data_index, 'mbsNum'] += mbsNum

                    if "Entropy" in metrics_list:
                        raw_data.loc[data_index, 'kmerAll'] += kmerAll
                        for i in range(self.K):
                            for j in range(self.K):
                                index = (i) * self.K + j
                                raw_data.loc[data_index, f'kmer{i}_{j}'] += kmer_list[index]
                except:
                    new_data = {'chrom': chrom, 'pos': pos, 'real': real, 'refer': refer, 'nReads': nReads, 'mBase': mBase,
                                'tBase': tBase, 'cBase': cBase, 'K4plus': K4plus, 'nDR': nDR, 'nMR': nMR}
                    if "MHL" in metrics_list:
                        for i in range(self.minK - 1, self.maxK):
                            new_data[f'methKmers{i}'] = methKmers_list[i]
                            new_data[f'totalKmers{i}'] = totalKmers_list[i]

                    if "MBS" in metrics_list:
                        new_data['mbsNum'] = mbsNum

                    if "Entropy" in metrics_list:
                        for i in range(self.K):
                            for j in range(self.K):
                                index = (i) * self.K + j
                                new_data[f'kmer{i}_{j}'] = kmer_list[index]
                                new_data['kmerAll'] = kmerAll

                    raw_data.loc[data_index] = new_data

        output_file = os.path.join(self.outputDir, self.tag + ".csv")
        raw_data.to_csv(output_file, index=False, header=True, sep='\t')

def main(args):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + "Run summaryBySNP start!")
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

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + "Run summaryBySNP end!")

if __name__ == '__main__':
    args = Namespace()
    args.epibedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRR1045636.epibed.gz"
    args.cpgPath = "/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg19/hg19_CpG.gz"
    args.fastaPath = "/sibcb2/bioinformatics/iGenome/Bismark/hg19/hg19.fa"
    # args.region = "chr1:10455-12640322"
    # args.bedPath = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/bed_file/test.bed"
    args.outputDir = "/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir"
    args.tag = "summaryBySNP.test"
    args.metrics = "MM PDR CHALM MHL MCR MBS Entropy"
    # args.metrics = "Entropy"
    args.minK = 1
    args.maxK = 10
    args.K = 2
    args.cutReads = True
    args.strand = "both"
    args.k4Plus = 2
    args.cpgCov = 10
    args.r2Cov = 10
    main(args)
