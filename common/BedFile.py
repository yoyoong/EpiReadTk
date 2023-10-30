import os, sys
# from common import Region
from common.Region import Region

class BedFile:
    def __init__(self, bedPath: str, base_type: int):
        self.bed_path = bedPath
        self.bed_name = os.path.splitext(os.path.basename(bedPath))[0]
        self.base_type = base_type

    def get_region_list(self):
        self.region_list = list()
        with open(self.bed_path) as file:
            for line in file.readlines():
                line = line.rstrip()
                if "\t" in line:
                    region_str = line.split("\t")
                elif " " in line:
                    region_str = line.split(" ")
                else:
                    print('Input bed file format error!')
                    sys.exit()

                if len(region_str) < 3:
                    print('Input bed file format error!')
                    sys.exit()

                if self.base_type == 0:
                    region = Region.init(region_str[0], region_str[1], region_str[2])
                elif self.base_type == 1:
                    region = Region.init(region_str[0], region_str[1] - 1, region_str[2])
                self.region_list.append(region)

        return self.region_list