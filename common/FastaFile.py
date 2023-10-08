from common import Region
import pysam


class FastaFile:
    def __init__(self, fastaPath: str):
        self.fasta_path = fastaPath

    def query_by_region(self, chrom: str, start: int, end: int):
        '''Query fasta file by region.
        Return:
            fasta_info: a string of base information.
        '''
        self.fasta_info = ""
        read_fa = pysam.FastaFile(self.fasta_path)
        self.fasta_info = read_fa.fetch(chrom, start, end)
        read_fa.close()
        return self.fasta_info