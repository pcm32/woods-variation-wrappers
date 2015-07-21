__author__ = 'pmoreno'

from TabixQuery.core import TabixQuery
import sys
import re

'''
Reads the tabix indexed vcf.bgzp from the first argument, and the pairs of Chromosome start from the
stdin
'''
if __name__ == '__main__':

    tabix_querier = TabixQuery(sys.argv[1])

    print "\t".join(["Chrom", "Pos", "Reference", "Alternate", "MAF_EuropeanAmerican"])
    for line in sys.stdin:
        chrom, position = line.split()
        for vcfEntry in tabix_querier.query(chrom,int(position),int(position)+1):
            match = re.search('MAF=([0\.]{0,2}\d+),([0\.]{0,2}\d+),([0\.]{0,2}\d+);', vcfEntry.info)
            if match is not None:
                print "\t".join([str(vcfEntry.chrom), str(vcfEntry.pos), vcfEntry.reference, vcfEntry.alt, str(match.group(1))])
