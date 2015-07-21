__author__ = 'pmoreno'

from pysam import TabixFile

class VCFEntry(object):

    def __init__(self, pysam_tabix_row):
        pysam_tabix_tuple = pysam_tabix_row.split('\t')
        self.chrom = pysam_tabix_tuple[0]
        self.pos = pysam_tabix_tuple[1]
        self.id = pysam_tabix_tuple[2]
        self.reference = pysam_tabix_tuple[3]
        self.alt = pysam_tabix_tuple[4]
        self.qual = pysam_tabix_tuple[5]
        self.filter = pysam_tabix_tuple[6]
        self.info = pysam_tabix_tuple[7]






class TabixQuery(object):

    def __init__(self, bgzipVCFPath):
        self.bzgipVCFPath = bgzipVCFPath
        self.reader = TabixFile(bgzipVCFPath)

    def query(self, chrom, start, stop):
        vcfEntries = []
        for row in self.reader.fetch(chrom,start,stop):
            vcfEntries.append(VCFEntry(row))
        return vcfEntries

    def __query(self, chrom, start):
        return self.query(chrom,start,start+1)