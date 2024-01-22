import argparse
import datetime
import gzip
import itertools
import math
import os
import sys
from pathlib import Path

import pandas as pd
from tqdm.autonotebook import tqdm

from main import HelpMessage
from common.BedFile import BedFile
from common.Constant import CODE, CPG_DICT, SNP_INSERT_DICT, SNP_REPLACE_ATCG_DICT
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile
from common.FastaFile import FastaFile
from common.Namespace import Namespace
from common.Region import Region


class SummaryBySNP:
    def __init__(self, args):
        self.epibedFile = EpibedFile(args.epibedPath)
        self.cpgFile = CpGFile(args.cpgPath)
        self.fastaFile = FastaFile(args.fastaPath)
        self.region = Region(args.region) if args.region is not None else None
        self.bedFile = BedFile(args.bedPath, 0) if args.bedPath is not None else None
        self.outputDir = args.outputDir
        self.tag = args.tag
        self.metrics_list = args.metrics.split(" ")
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
        # assert 3 <= self.K <= 5, 'K: the default is 4 and values must be between 3 and 5'
        # assert self.region or self.bedFile, 'you should input bedPath or region'
        assert (not self.region) or (not self.bedFile), 'you should only input bedPath or region'
        assert self.strand == 'both' or self.strand == 'plus' or self.strand == 'minus', 'strand should be both, plus or minus'

    def get_columns_list(self):
        column_list = ["chrom", "pos", "refer", "real", "nReads", "mBase", "cBase", "tBase", "K4plus", "nDR", "nMR"]
        if "MHL" in self.metrics_list:
            for i in range(self.minK - 1, self.maxK):
                column_list.append(f'methKmers{i}')
                column_list.append(f'totalKmers{i}')
        if "MBS" in self.metrics_list:
            column_list.append('mbsNum')
        if "Entropy" in self.metrics_list:
            for i in range(args.K):
                for j in range(args.K):
                    column_list.append(f'kmer{i}_{j}')
            column_list.append('kmerAll')
        return column_list

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

    def process_line(self, info:list, raw_data:pd.DataFrame):
        if len(info) < 9:
            return
        start_pos = int(info[1]) + 1  # 0-base to 1-base
        read_strand = info[5]
        if self.strand != 'both':
            if self.strand != read_strand:
                return

        # get the snp number in this info
        snp_num = 0
        for s in info[8]:
            if s in SNP_REPLACE_ATCG_DICT.keys():  # 仅考虑ATCG replacement
                snp_num += 1
        if snp_num == 0:
            return

        cpg_info = self.separate_rle_string(info[6])
        snp_info = self.separate_rle_string(info[8])

        # get the snp position and label list
        snp_pos_list = []
        snp_label_list = []
        snp_move_length = 0
        snp_insert_length = 0
        for i in range(len(snp_info)):
            if snp_info[i].isalpha():
                if snp_info[i] in SNP_REPLACE_ATCG_DICT.keys():  # 仅考虑ATCG replacement
                    if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                        for s in range(int(snp_info[i + 1])):
                            position = start_pos + snp_move_length - snp_insert_length + s
                            snp_pos_list.append(position)
                            snp_label_list.append(snp_info[i])
                        snp_move_length += int(snp_info[i + 1])
                    else:
                        position = start_pos + snp_move_length - snp_insert_length
                        snp_pos_list.append(position)
                        snp_label_list.append(snp_info[i])
                        snp_move_length += 1
                elif snp_info[i] in SNP_INSERT_DICT.keys():
                    if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                        snp_insert_length += int(snp_info[i + 1])
                    else:
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
            return

        # get the cpg label string, '1'-methylated, '0'-unmethylated, '-' is not cover
        chrom = f"chr{info[0]}" if not str(info[0]).startswith("chr") else str(info[0])
        region = Region.init(chrom, cpg_pos_list[0], cpg_pos_list[-1] + 1)
        try:
            cpg_full_pos_list = self.cpgFile.query_by_region(region)
        except:
            print(f"{region.format_string} are not exist in cpg reference file, continue...")
            return
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
            pos = snp_pos
            real = snp_label_list[i]
            refer = self.fastaFile.query_by_region(info[0], snp_pos, snp_pos + 1).upper()
            nReads = 1  # 总read个数
            mBase = cpg_label_list.count('1')  # 甲基化位点个数
            tBase = len(cpg_raw_label_list)  # 总位点个数
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
            if "MHL" in self.metrics_list:
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
            if "MBS" in self.metrics_list:
                if len(cpg_label_list) >= self.K:
                    cpg_label_split_null = cpg_label_str.split('-')
                    for item in cpg_label_split_null:
                        cpg_label_split = list(filter(None, item.split('0')))
                        mbsNum = sum([len(x) ** 2 for x in cpg_label_split]) / len(cpg_label_list) ** 2

            kmer_list = [0] * (self.K * self.K)  # every kmer pattern number
            kmerAll = 0
            if "Entropy" in self.metrics_list:
                if len(cpg_label_list) >= self.K:
                    for i in range(len(cpg_label_str) - self.K + 1):
                        kmer = cpg_label_str[i: i + self.K]
                        if '-' not in kmer:
                            kmer_list[int(kmer, 2)] += 1
                            kmerAll += 1

            data_index = chrom + str(pos) + real + refer
            try:  # try success mean data_index exist
                raw_data.loc[data_index, 'nReads'] += nReads
                raw_data.loc[data_index, 'mBase'] += mBase
                raw_data.loc[data_index, 'tBase'] += tBase
                raw_data.loc[data_index, 'cBase'] += cBase
                raw_data.loc[data_index, 'K4plus'] += K4plus
                raw_data.loc[data_index, 'nDR'] += nDR
                raw_data.loc[data_index, 'nMR'] += nMR
                if "MHL" in self.metrics_list:
                    for i in range(self.minK - 1, self.maxK):
                        raw_data.loc[data_index, f'methKmers{i}'] += methKmers_list[i]
                        raw_data.loc[data_index, f'totalKmers{i}'] += totalKmers_list[i]

                if "MBS" in self.metrics_list:
                    raw_data.loc[data_index, 'mbsNum'] += mbsNum

                if "Entropy" in self.metrics_list:
                    raw_data.loc[data_index, 'kmerAll'] += kmerAll
                    for i in range(self.K):
                        for j in range(self.K):
                            index = (i) * self.K + j
                            raw_data.loc[data_index, f'kmer{i}_{j}'] += kmer_list[index]
            except:
                new_data = {'chrom': chrom, 'pos': pos, 'real': real, 'refer': refer, 'nReads': nReads,
                            'mBase': mBase, 'tBase': tBase, 'cBase': cBase, 'K4plus': K4plus, 'nDR': nDR,
                            'nMR': nMR}
                if "MHL" in self.metrics_list:
                    for i in range(self.minK - 1, self.maxK):
                        new_data[f'methKmers{i}'] = methKmers_list[i]
                        new_data[f'totalKmers{i}'] = totalKmers_list[i]

                if "MBS" in self.metrics_list:
                    new_data['mbsNum'] = mbsNum

                if "Entropy" in self.metrics_list:
                    for i in range(self.K):
                        for j in range(self.K):
                            index = (i) * self.K + j
                            new_data[f'kmer{i}_{j}'] = kmer_list[index]
                            new_data['kmerAll'] = kmerAll

                raw_data.loc[data_index] = new_data

    def summary_by_region(self, region: Region):
        self.epibed_info = self.epibedFile.query_by_region(region)
        column_list = self.get_columns_list()
        raw_data = pd.DataFrame(columns=column_list)

        for info in self.epibed_info:
            self.process_line(info, raw_data)

        return raw_data

    def calculate_all(self):
        column_list = self.get_columns_list()
        raw_data = pd.DataFrame(columns=column_list)

        count = 0
        agg_dict = {} # 合并raw_data时每一列的处理规则
        for column in raw_data.columns:
            if column in ["chrom", "pos", "refer", "real"]:
                agg_dict[column] = 'first'
            else:
                agg_dict[column] = 'sum'

        first_flag = True
        for line in tqdm(gzip.open(self.epibedFile.epibed_path, 'rb'), "Read EpiRead"):
            count += 1
            if count % 100000 == 0: # 每10万行重置raw_data，避免raw_data过大
                if first_flag:
                    self.write_to_file(raw_data, mode="w")
                    first_flag = False
                else:
                    self.write_to_file(raw_data, mode="a")
                raw_data = pd.DataFrame(columns=column_list)

            info = line.decode().split('\t')
            # if info[0] != '3':
            #     continue
            # else:
            #     if info[1] < 90449414:
            #         continue
            self.process_line(info, raw_data)

    def write_to_file(self, raw_data: pd.DataFrame, mode: str = 'w') -> None:
        head_base_list = ["chrom", "pos", "refer", "real", "nReads", "mBase", "cBase", "tBase", "K4plus", "nDR", "nMR"]
        output_column_list = head_base_list + self.metrics_list
        output_filename = os.path.join(self.outputDir, f"{self.tag}.csv")

        # filter out only have one snp position
        raw_data['duplicated'] = raw_data.duplicated(['pos'], keep=False)
        output_data = raw_data[raw_data['duplicated'] == True].sort_values(by=["chrom", "pos", "refer", "real"])
        for index in output_data.index:
            if "MM" in self.metrics_list:
                if output_data.loc[index, "nReads"] < self.cpgCov:
                    output_data.loc[index, "MM"] = float('nan')
                else:
                    MM = round((output_data.loc[index, "mBase"] / output_data.loc[index, "tBase"]), 4)
                    output_data.loc[index, "MM"] = MM
            if "PDR" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "PDR"] = float('nan')
                else:
                    PDR = round(output_data.loc[index, "nDR"] / output_data.loc[index, "K4plus"], 4)
                    output_data.loc[index, "PDR"] = PDR
            if "CHALM" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "CHALM"] = float('nan')
                else:
                    CHALM = round(output_data.loc[index, "nMR"] / output_data.loc[index, "K4plus"], 4)
                    output_data.loc[index, "CHALM"] = CHALM
            if "MHL" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "MHL"] = float('nan')
                else:
                    temp = 0.0
                    weight = 0
                    for j in range(self.minK - 1, self.maxK):
                        methKmers = output_data.loc[index, f"methKmers{j}"]
                        totalKmers = output_data.loc[index, f"totalKmers{j}"]
                        if totalKmers < 1:
                            continue
                        temp += methKmers / totalKmers * (j + 1);
                        weight += j + 1
                    MHL = round(temp / weight)
                    output_data.loc[index, "MHL"] = MHL
            if "MCR" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "MCR"] = float('nan')
                else:
                    MCR = round(output_data.loc[index, "cBase"] / output_data.loc[index, "tBase"], 4)
                    output_data.loc[index, "MCR"] = MCR
            if "MBS" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "MBS"] = float('nan')
                else:
                    MBS = round(output_data.loc[index, "mbsNum"] / output_data.loc[index, "K4plus"], 4)
                    output_data.loc[index, "MBS"] = MBS
            if "Entropy" in self.metrics_list:
                if output_data.loc[index, "K4plus"] < self.k4Plus:
                    output_data.loc[index, "Entropy"] = float('nan')
                else:
                    temp = 0.0
                    for m in range(self.K):
                        for n in range(self.K):
                            kmer_num = output_data.loc[index, f"kmer{m}_{n}"]
                            kmerAll_num = output_data.loc[index, f"kmerAll"]
                            if kmer_num < 1 or kmerAll_num < 1:
                                continue
                            temp += kmer_num / kmerAll_num * math.log(kmer_num / kmerAll_num, 2)
                    Entropy = round(- 1 / self.K * temp, 4)
                    output_data.loc[index, "Entropy"] = Entropy

        if output_data.shape[0] == 0:
            return
        if mode == 'w':
            output_data.to_csv(output_filename, columns=output_column_list, index=False, header=True, sep='\t', mode=mode)
        else:
            output_data.to_csv(output_filename, columns=output_column_list, index=False, header=False, sep='\t', mode=mode)

def main(args):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f ') + "Run summaryBySNP start!")
    summaryBySNP = SummaryBySNP(args)
    summaryBySNP.check_args()

    if summaryBySNP.region is not None:
        raw_data = summaryBySNP.summary_by_region(summaryBySNP.region)
        summaryBySNP.write_to_file(raw_data, mode="w")
    elif summaryBySNP.bedFile is not None:
        region_list = summaryBySNP.bedFile.get_region_list()
        first_flag = True
        for region in tqdm(region_list, "Summary by region"):
            raw_data = summaryBySNP.summary_by_region(region)
            if first_flag:
                summaryBySNP.write_to_file(raw_data, mode="w")
                first_flag = False
            else:
                summaryBySNP.write_to_file(raw_data, mode="a")
    else:
        summaryBySNP.calculate_all()

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f ') + "Run summaryBySNP end!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='summaryBySNP')
    help = HelpMessage()
    parser.add_argument('--epibedPath', type=str, required=False, help=help.epibedPath,
                        default="/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/Z000000PX.epibed.gz") # require
    parser.add_argument('--cpgPath', type=str, required=False, help=help.cpgPath,
                        default="/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg38/hg38_CpG.gz") # require
    parser.add_argument('--fastaPath', type=str, required=False, help=help.fastaPath,
                        default="/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/0.refer_fastq/hg38.fa") # require
    parser.add_argument('--region', type=str, required=False, help=help.region)
    parser.add_argument('--bedPath', type=str, required=False, help=help.bedPath)
    parser.add_argument('--outputDir', type=str, required=False, help=help.outputDir,
                        default="/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir") # require
    parser.add_argument('--tag', type=str, required=False, help=help.tag, default="Z000000PX_18") # require
    parser.add_argument('--metrics', type=str, required=False, help=help.metrics,
                        default="MM PDR CHALM MHL MCR MBS Entropy")
    parser.add_argument('--minK', type=int, required=False, help=help.minK, default='1')
    parser.add_argument('--maxK', type=int, required=False, help=help.maxK, default='10')
    parser.add_argument('--K', type=int, required=False, help=help.K, default='4')
    parser.add_argument('--cutReads', required=False, action='store_true', help=help.cutReads)
    parser.add_argument('--strand', type=str, required=False, help=help.strand, default='both')
    parser.add_argument('--k4Plus', type=int, required=False, help=help.k4Plus, default='5')
    parser.add_argument('--cpgCov', type=int, required=False, help=help.cpgCov, default='5')
    parser.add_argument('--r2Cov', type=int, required=False, help=help.r2Cov, default='5')
    args = parser.parse_args()

    main(args)