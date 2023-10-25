
class Region:
    def __init__(self, region: str):
        self.chr = region.split(':')[0]
        self.start = int(region.split(':')[1].split('-')[0])
        self.end = int(region.split(':')[1].split('-')[1])

    @classmethod
    def init(self, chrom: str, start: int, end: int):
        region_str = chrom + ":" + str(start) + "-" + str(end)
        return self(region_str)

    @property
    def len(self) -> int:
        return self.end - self.start

    @property
    def format_string(self) -> str:
        return self.chr + ":" + str(self.start) + "-" + str(self.end)

    @property
    def filename_string(self) -> str:
        return self.chr + "_" + str(self.start) + "_" + str(self.end)