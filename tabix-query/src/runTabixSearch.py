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

    print "\t".join(["Chrom", "Pos", "Reference", "Alternate", "MAF_EuropeanAmerican", "EA_Ref_Count", "EA_Alt_Count"])
    for line in sys.stdin:
        chrom, position = line.split()
        for vcfEntry in tabix_querier.query(chrom,int(position)-1,int(position)):
            #match = re.search('MAF=([0\.]{0,2}\d+),([0\.]{0,2}\d+),([0\.]{0,2}\d+);', vcfEntry.info)
            match = re.search('MAF=(\d+\.\d+),(\d+\.\d+),(\d+\.\d+);', vcfEntry.info)
            maf = None
            if match is not None:
                maf = match.group(1)
            match = re.search('EA_AC=(\d+),(\d+);', vcfEntry.info)
            alt_count = None
            ref_count = None
            if match is not None:
                alt_count = match.group(1)
                ref_count = match.group(2)
            print "\t".join([str(vcfEntry.chrom), str(vcfEntry.pos), vcfEntry.reference, vcfEntry.alt,
                             str(maf), str(ref_count), str(alt_count)])
