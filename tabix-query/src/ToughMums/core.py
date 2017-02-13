__author__ = 'pmoreno'


class ToughMumsEntry(object):

    def __init__(self, entryLine):
        self.completeLine = entryLine
        self.Chromosome, self.Position, self.Change, self.Cohort_Allele_Count, \
        self.Cohort_Allele_Frequency, self.T1000G_Allele_Count, \
        self.T1000G_Allele_Frequency, self.Chi_Squared_Statistic, \
        self.P_Value, self.Bonferroni_Correction_005_Significance, \
        self.Gene, self.Effect, self.cDNA_pos, self.codon_pos, \
        self.Protein_pos, self.AA_change, self.Grantham_score, \
        self.dbsnp, self.PolyPhen, self.SIFT, self.Protein_Domain, \
        self.Clinic_Sig, self.Canonical_Trans, self.GERP, self.PHYLOP10 \
            = entryLine.split(sep='\t')


class ToughMumsResultReader(object):

    def __init__(self, pathToResultFile):
        self.reader = open(pathToResultFile, 'r')
        self.reader.readline()

    def next(self):
        return ToughMumsEntry(self.reader.readline().rstrip('\n'))

    def close(self):
        self.reader.close()


class ToughMumsWriter(object):

    def __init__(self, path_out, header):
        self.writer = open(path_out, 'w')
        self.writer.write(header+'\n')

    def write(self, toughMumsEntry):
        self.writer
