__author__ = 'pmoreno'

from TabixQuery.core import TabixQuery
import sys
import fileinput
import re

'''
Reads the tabix indexed vcf.bgzp from the first argument, and the pairs of Chromosome start from the
stdin
'''
if __name__ == '__main__':

    tabix_querier = TabixQuery(sys.argv[1])

    for line in fileinput.input():
        chrom, position = line.split(sep=" ")
        for vcfEntry in tabix_querier.query(chrom,position):
            match = re.search('MAF=([0\.]{0,2}\d+),([0\.]{0,2}\d+),([0\.]{0,2}\d+);', vcfEntry.info)
            print "\t".join([vcfEntry.chrom, vcfEntry.pos, vcfEntry.reference, vcfEntry.alt, match.group(1)])+"\n"