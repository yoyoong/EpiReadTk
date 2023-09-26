
class Region:
    def __init__(self, region: str, base_type: int):
        self.chr = region.split(':')[0]
        self.base_type = base_type
        if self.base_type == 0:
            self.start = int(region.split(':')[1].split('-')[0]) - 1
        elif self.base_type == 1:
            self.start = int(region.split(':')[1].split('-')[0])
        self.end = int(region.split(':')[1].split('-')[1])

    @property
    def len(self) -> int:
        return self.end - self.start

    @property
    def format_string(self) -> str:
        return self.chr + ":" + str(self.start) + "-" + str(self.end)

    @property
    def filename_string(self) -> str:
        return self.chr + "_" + str(self.start) + "_" + str(self.end)