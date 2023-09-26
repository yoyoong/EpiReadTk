import tabix
from common import Region
import itertools
import numpy as np

CPG_LABEL = ['M', 'U']
NOMESEQ_CPG_LABEL = ['O', 'S']
SNP_LABEL = ['D', 'd', 'A', 'T', 'C', 'G', 'R', 'Y', 'a', 't', 'c', 'g', 'i']

class EpibedFile:
    def __init__(self, epibedPath: str):
        self.epibed_path = epibedPath
        self.epibed_tabix = tabix.open(epibedPath)

    def query_by_region(self, region: Region):
        '''Query epibed file to by region.
        Return:
            epibed_info: a dictionary which key is read name and value is a list include 1 or 2(usually) epibed line.
        '''
        self.epibed_info = dict()
        records = self.epibed_tabix.query(region.chr, region.start, region.end)
        for item in records:
            if len(item) < 9:
                continue
            item[1] = int(item[1]) + 1  # 0-base to 1-base
            read_name = item[3]
            epibed_read_info = self.epibed_info.get(read_name)
            if epibed_read_info is None:
                epibed_read_info = []
            epibed_read_info.append(item)
            self.epibed_info[read_name] = epibed_read_info

        return self.epibed_info

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

    def get_cpg_snp_position(self):
        '''Get the genome position of CpG and SNP site.

        Return:
            cpg_snp_position: a list of genome position of CpG and SNP site.
        '''
        self.cpg_snp_position = set()
        for read in self.epibed_info.values():
            for line in read:
                start_pos = int(line[1])
                strand = line[5]
                cpg_info = self.separate_rle_string(line[6])
                snp_info = self.separate_rle_string(line[8])

                cpg_move_length = 0
                for i in range(len(cpg_info)):
                    if cpg_info[i].isalpha():
                        if cpg_info[i] in CPG_LABEL:
                            if strand == "+": # reads from OT/CTOT(+) strands, methylation site is in C→T substitution
                                self.cpg_snp_position.add(start_pos + cpg_move_length)
                            else: # reads from OB/CTOB (-) strands, methylation site is in G→A substitution
                                self.cpg_snp_position.add(start_pos + cpg_move_length - 1)

                        if i + 1 < len(cpg_info) and cpg_info[i + 1].isdigit():
                            cpg_move_length += int(cpg_info[i + 1])
                        else:
                            cpg_move_length += 1

                snp_move_length = 0
                for i in range(len(snp_info)):
                    if snp_info[i].isalpha():
                        if snp_info[i] in SNP_LABEL:
                            self.cpg_snp_position.add(start_pos + snp_move_length)

                        if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                            snp_move_length += int(snp_info[i + 1])
                        else:
                            snp_move_length += 1
                            # if strand == "+":
                            #     self.cpg_snp_position.add(start_pos if i == 0 else start_pos + snp_move_length)
                            # else:
                            #     self.cpg_snp_position.add(start_pos if i == 0 else start_pos + snp_move_length - 1)

        self.cpg_snp_position = list(sorted(self.cpg_snp_position))
        return self.cpg_snp_position

    def build_cpg_snp_matrix(self):
        '''Build a matrix include both CpG and SNP information.

        Return:
            cpg_snp_matrix: a matrix which row index is position and every row is a read's CpG and SNP information.In
            epiBED format, CpG information locate in columns 7 and SNP information in columns 9. The matrix replaces the
            labels in the CpG and SNP information with numbers, the meanings of the different numbers are as follows:
            0: ignore site
            -1: Base in reference was deleted in read at that location(D/d)
            10: Unmethylated CpG(U)
            11: Methylated CpG(M)
            21/22/23/24/25/26: SNP base seen in read relative to reference(A/T/C/G/R/Y); R and Y represent A/G and C/T
            31/32/33/34/35: Inserted base included in read, but no in reference(a/t/c/g/i); i is used as a placeholder
        '''
        pos_num = len(self.cpg_snp_position)
        read_num = len(self.epibed_info)
        self.cpg_snp_matrix = np.zeros((read_num, pos_num), dtype='int')

        row = 0
        for read in self.epibed_info.values():
            for line in read:
                start_pos = int(line[1])
                strand = line[5]
                cpg_info = self.separate_rle_string(line[6])
                snp_info = self.separate_rle_string(line[8])

                cpg_move_length = 0
                for i in range(len(cpg_info)):
                    if cpg_info[i].isalpha():
                        if cpg_info[i] in CPG_LABEL:
                            if strand == "+":  # reads from OT/CTOT(+) strands, methylation site is in C→T substitution
                                position = start_pos + cpg_move_length
                            else:  # reads from OB/CTOB (-) strands, methylation site is in G→A substitution
                                position = start_pos + cpg_move_length - 1
                            col = self.cpg_snp_position.index(position)

                            if cpg_info[i] == "U":
                                self.cpg_snp_matrix[row][col] = 10
                            elif cpg_info[i] == "M":
                                self.cpg_snp_matrix[row][col] = 11

                        if i + 1 < len(cpg_info) and cpg_info[i + 1].isdigit():
                            cpg_move_length += int(cpg_info[i + 1])
                        else:
                            cpg_move_length += 1

                snp_move_length = 0
                for i in range(len(snp_info)):
                    if snp_info[i].isalpha():
                        if snp_info[i] in SNP_LABEL:
                            position = start_pos + snp_move_length
                            col = self.cpg_snp_position.index(position)

                            if snp_info[i] == "D" or snp_info[i] == "d":
                                self.cpg_snp_matrix[row][col] = -1
                            elif snp_info[i] == "A":
                                self.cpg_snp_matrix[row][col] = 21
                            elif snp_info[i] == "T":
                                self.cpg_snp_matrix[row][col] = 22
                            elif snp_info[i] == "C":
                                self.cpg_snp_matrix[row][col] = 23
                            elif snp_info[i] == "G":
                                self.cpg_snp_matrix[row][col] = 24
                            elif snp_info[i] == "R":
                                self.cpg_snp_matrix[row][col] = 25
                            elif snp_info[i] == "Y":
                                self.cpg_snp_matrix[row][col] = 26
                            elif snp_info[i] == "a":
                                self.cpg_snp_matrix[row][col] = 31
                            elif snp_info[i] == "t":
                                self.cpg_snp_matrix[row][col] = 32
                            elif snp_info[i] == "c":
                                self.cpg_snp_matrix[row][col] = 33
                            elif snp_info[i] == "g":
                                self.cpg_snp_matrix[row][col] = 34
                            elif snp_info[i] == "i":
                                self.cpg_snp_matrix[row][col] = 35

                        if i + 1 < len(snp_info) and snp_info[i + 1].isdigit():
                            snp_move_length += int(snp_info[i + 1])
                        else:
                            snp_move_length += 1
            row += 1

        return self.cpg_snp_matrix