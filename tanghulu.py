import os
from pathlib import Path
from common.Dict import Dict
from common.Region import Region
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile
from common.Constant import CODE, CPG_DICT, SNP_DICT
import matplotlib.pyplot as plt
import numpy as np


SHIFT = 500

class Tanghulu:
    def __init__(self, args):
        self.epibedFile = EpibedFile(args.epibedFile)
        self.cpgFile = CpGFile(args.cpgFile, SHIFT)
        self.region = Region(args.region, 0)
        self.outputDir = args.outputDir
        self.tag = args.tag
        self.outFormat = args.outFormat
        self.outcut = args.outcut

    def check_args(self):
        assert Path(os.path.dirname(self.outputDir)).exists(), 'output directory does not exist'
        # assert self.outcut > self.region.len, f"The region is larger than {args.outcut}"

    def build_epibed_matrix(self):
        self.epibed_info = self.epibedFile.query_by_region(self.region)
        self.cpg_snp_position = self.epibedFile.get_cpg_snp_position()
        self.cpg_snp_matrix, self.strand_list = self.epibedFile.build_cpg_snp_matrix()

    def paint_tanghulu_plot(self):
        cpg_pos_list = self.cpgFile.query_by_region(self.region)
        ref_pos = self.cpg_snp_position
        show_pos = np.arange(len(ref_pos))
        pos_num = len(ref_pos)
        read_num = len(self.epibed_info)

        # figure size
        plt.clf()
        fig = plt.figure()
        fig.set_size_inches(pos_num + 3, read_num)

        # x-axis
        ax = plt.axes([0.1, 0.1, 0.7, 0.8])
        ax.plot(show_pos, [0] * pos_num, c="#c0c0c0", zorder=1)
        for pos in range(pos_num):
            if ref_pos[pos] in cpg_pos_list: # draw grey circle in cpg site
                ax.scatter(pos, [0], c="#c0c0c0", s=100, zorder=2)
            ax.text(pos, -0.5, ref_pos[pos], rotation=45, ha='center', rotation_mode='anchor', color="grey", zorder=3)
        plt.xticks([])
        plt.xlabel(xlabel='Genomic position', fontdict={'fontsize': 15})

        # y-axis
        plt.yticks(range(read_num + 1))

        # tanghulu
        y_pos = 1
        for i in range(read_num):
            # tanghulu color
            strand = self.strand_list[i]
            color = 'black' if strand == "+" else 'blue'

            # draw line
            status_array = self.cpg_snp_matrix[i]
            if np.all(status_array == 0):
                continue
            first_index = np.nonzero(np.array(status_array))[0].tolist()[0] # the first non-zero index
            last_index = np.nonzero(np.array(status_array))[0].tolist()[-1] # the last non-zero index
            line_start_x = show_pos[first_index]
            line_end_x = show_pos[last_index]
            ax.plot([line_start_x, line_end_x], [y_pos, y_pos], c=color, zorder=1)

            for j in range(pos_num):
                status = status_array[j]

                # draw white circle as placehold
                cpg_index = np.in1d(status_array, [value for value in CPG_DICT.values()])
                cpg_placehold_circle = show_pos[cpg_index]
                ax.scatter(cpg_placehold_circle, [y_pos] * len(cpg_placehold_circle), c='white', s=100, zorder=1)

                # draw cpg information
                if status in CPG_DICT.values():
                    ax.scatter([show_pos[j]], [y_pos], c='white', edgecolors=color, s=100, zorder=2)
                    if status == CODE.METHYLATED.value:
                        ax.scatter([show_pos[j]], [y_pos], c=color, edgecolors=color, s=100, zorder=3)

                # draw white circle as placehold
                snp_index = np.in1d(status_array, [value for value in SNP_DICT.values()])
                snp_placehold_circle = show_pos[snp_index]
                ax.scatter(snp_placehold_circle, [y_pos] * len(snp_placehold_circle), c='white', s=100, zorder=1)

                # draw snp information
                if status in SNP_DICT.values():
                    marker = ''
                    if status == CODE.DELETION.value:
                        marker = CODE.DELETION.label1
                    elif status == CODE.A_REPLACE.value:
                        marker = CODE.A_REPLACE.label
                    elif status == CODE.T_REPLACE.value:
                        marker = CODE.T_REPLACE.label
                    elif status == CODE.C_REPLACE.value:
                        marker = CODE.C_REPLACE.label
                    elif status == CODE.G_REPLACE.value:
                        marker = CODE.G_REPLACE.label
                    elif status == CODE.AorG_REPLACE.value:
                        marker = CODE.AorG_REPLACE.label
                    elif status == CODE.CorT_REPLACE.value:
                        marker = CODE.CorT_REPLACE.label
                    elif status == CODE.A_INSERT.value:
                        marker = CODE.A_INSERT.label
                    elif status == CODE.T_INSERT.value:
                        marker = CODE.T_INSERT.label
                    elif status == CODE.C_INSERT.value:
                        marker = CODE.C_INSERT.label
                    elif status == CODE.G_INSERT.value:
                        marker = CODE.G_INSERT.label
                    elif status == CODE.UNKNOWN_INSERT.value:
                        marker = CODE.UNKNOWN_INSERT.label

                    ax.scatter([show_pos[j]], [y_pos], c='red', marker=f"${marker}$", s=100, zorder=4)

            y_pos += 1

        # title
        ax.set_title(f'{self.region.format_string} ({self.epibedFile.epibed_name})', fontdict = {'fontsize': 20})

        # legend
        ax.axis(ymin=(-1))
        plt.scatter([], [], c='black', marker='_', s=50, label='OT/CTOT(+) strand')
        plt.scatter([], [], c='blue', marker='_', s=50, label='OB/CTOB(-) strand')
        plt.scatter([], [], c='black', s=50, label='Methylated')
        plt.scatter([], [], c='white', edgecolors='black', s=50, label='Unmethylated')
        plt.scatter([], [], c='red', marker='$\mathtt{D/d}$', s=150, label='Deleted base')
        plt.scatter([], [], c='red', marker='$\mathtt{A/T/C/G/R/Y}$', s=1600, label='SNP base')
        plt.scatter([], [], c='red', marker='$\mathtt{a/t/c/g/i}$', s=1200, label='Inserted base')
        plt.legend(prop={"size": 11}, loc=7, bbox_to_anchor=(1.22, 0.5), markerscale = 2, handlelength=7, handleheight=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        output_fig = os.path.join(self.outputDir, self.tag + "." + self.outFormat)
        fig.savefig(output_fig)
        # plt.close()
        plt.show()

def main(args):

    tanghulu = Tanghulu(args)
    tanghulu.check_args()
    tanghulu.build_epibed_matrix()
    tanghulu.paint_tanghulu_plot()

if __name__ == '__main__':
    args = Dict()
    args["epibedFile"] = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args["cpgFile"] = "/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg19/hg19_CpG.gz"
    args["region"] = "chr1:190027767-190029558"
    args["outputDir"] = "/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir"
    args["tag"] = "tanghulu.test"
    args["outFormat"] = "png"
    args["outcut"] = 2000
    main(args)
