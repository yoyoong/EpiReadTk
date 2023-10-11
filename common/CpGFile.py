import tabix
from common import Region

class CpGFile:
    def __init__(self, cpgPath: str, shift = 0):
        self.cpg_path = cpgPath
        self.shift = shift
        self.cpg_tabix = tabix.open(cpgPath)

    def query_by_region(self, region: Region):
        cpg_pos_list = []
        records = self.cpg_tabix.query(region.chr, region.start - self.shift, region.end + self.shift)
        for item in records:
            if len(item) < 3:
                continue
            cpg_pos_list.append(int(item[1]))
        return cpg_pos_list