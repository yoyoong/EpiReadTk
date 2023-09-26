import os
from pathlib import Path
from common.Dict import Dict
from common.Region import Region
from common.CpGFile import CpGFile
from common.EpibedFile import EpibedFile

SHIFT = 500

class Tanghulu:
    def __init__(self, args):
        self.epibedPath = EpibedFile(args.epibedPath)
        self.cpgPath = CpGFile(args.cpgPath, SHIFT)
        self.region = Region(args.region, 0)
        self.outputDir = args.outputDir
        self.tag = args.tag
        self.outFormat = args.outFormat
        self.outcut = args.outcut

    def check_args(self):
        assert Path(os.path.dirname(self.outputDir)).exists(), 'output directory does not exist'
        # assert self.outcut > self.region.len, f"The region is larger than {args.outcut}"

    def build_epibed_matrix(self):
        self.epibedPath.query_by_region(self.region)
        self.epibedPath.get_cpg_snp_position()
        self.epibedPath.build_cpg_snp_matrix()

    def paint_tanghulu(self):
        cpg_pos_list = self.cpgPath.query_by_region(self.region)
        print(cpg_pos_list)

def main(args):

    tanghulu = Tanghulu(args)
    tanghulu.check_args()
    tanghulu.build_epibed_matrix()

if __name__ == '__main__':
    args = Dict()
    args["epibedPath"] = "/sibcb2/bioinformatics2/hongyuyang/project/EpiReadTk/data/6.epibed/SRX1631736.epibed.gz"
    args["cpgPath"] = "/sibcb2/bioinformatics2/zhangzhiqiang/genome/CpG/hg19/hg19_CpG.gz"
    args["region"] = "chr1:546494-10547092"
    args["outputDir"] = "/sibcb2/bioinformatics2/hongyuyang/code/EpiReadTk/outputDir"
    args["tag"] = "tanghulu.test"
    args["outFormat"] = "pdf"
    args["outcut"] = 2000
    main(args)
