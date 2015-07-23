suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

options_list<-list(
  make_option(c("--cohortCounts","-c"),help="cohorts count file",),
  make_option(c("--output","-o"),help="ouput full path name"),
  make_option(c("--tabixResult","-t"),help="full path to tabix processed result"),
  make_option(c("--femaleCount","-f"),help="Number of females"),
  make_option(c("--maleCount","-m"),help="Number of males"),
  make_option(c("--ref1000GPath","-r"),help="Path to 1000 Genomes reference file")
  make_option(c("--bamUnseqResult","-u"),help="Path to file with the file names containing the positions that are not sequenced per sample"),
)

opt<-parse_args(OptionParser(option_list = options_list),positional_arguments = FALSE)

# we index all tables by a new column: chr-position-allele
# cohort columns: Chromosome      Position        Observed_Allele Allele_Count    Assumed_Total_Alleles
fread(opt$cohortCounts)->cohortCounts
cohortCounts[,Assumed_Total_Alleles:=opt$femaleCount+opt$maleCount,]
cohortCounts[Chromosome=="Y",Assumed_Total_Alleles:=opt$maleCount,]
cohortCounts[,Other_Alleles_Count:=Assumed_Total_Alleles-Allele_Count,]
cohortCounts[,AlleleKey:=paste(Chromosome,Position,Observed_Allele,sep="_"),]
setkey(cohortCounts,AlleleKey)

runChiSqrd<-function(cohortAllele,cohortOther,referenceAllele,referenceOther) {
  if(is.na(referenceAllele)) {
    return(c(NA_real_,NA_real_))
  }
  Xsq<-chisq.test(as.table(rbind(c(cohortAllele,cohortOther),c(referenceAllele,referenceOther))))
  return(c(Xsq$p.value,Xsq$statistic))
}

#' if a 1000 genomes path is given, then we merge our cohort with that data.
#' This should the the formatted 1000 genomes provided by Katie.
if("ref1000GPath" %in% names(opt)) {
    # 1000 Genomes formatted file columns
    # Chromosome      Position        Observed_Allele Allele_Count    Total_Alleles   Gene    Ref_Allele      Effect  cDNA pos        codon pos       Protein pos     AA change       Grantham score  dbsnp   PolyPhen        SIFT    Protein Domain  Clinic Sig      Canonical Trans GERP    PHYLOP100
    fread(opt$ref1000GPath)->reference
    reference[,AlleleKey:=paste(Chromosome,Position,Observed_Allele,sep="_"),]
    #setkey(reference,AlleleKey)
    setkeyv(reference,c('AlleleKey','Effect'))
    
    reference[,list(
      Obs_Allele_1000G=unique(Observed_Allele)[1],
      Allele_Count_1000G=unique(Allele_Count)[1],
      Total_Alleles_1000G=unique(Total_Alleles)[1],
      Gene=unique(Gene)[1],Ref_Allele_1000G=unique(Ref_Allele)[1],
      cDNA_pos=unique(get('cDNA pos'))[1],
      Codon_pos=unique(get('codon pos'))[1],
      Protein_pos=unique(get('Protein pos'))[1],
      AA_change=unique(get('AA change'))[1],
      Grantham_Score=unique(get('Grantham score'))[1],
      dbsnp=paste(unique(unlist(strsplit(dbsnp,split = ";"))),collapse = ";"), 
      PolyPhen=unique(PolyPhen)[1],SIFT=unique(SIFT)[1],
      Protein_Domain=paste(unique(unlist(strsplit(get('Protein Domain'),split='&'))),collapse = "&"),
      Clinic_Sig=paste(unique(get('Clinic Sig')),collapse="-"),
      Canonical_Trans=gsub(pattern = '(^-|-$)', replacement = "", x = paste(unique(get('Canonical Trans')),collapse="-")),
      GERP=unique(GERP)[1],PHYLOP100=unique(PHYLOP100)[1]),
    by=c('AlleleKey','Effect')]->ref_short
    
    setkey(ref_short,AlleleKey)
    
    ref_short[cohortCounts,allow.cartesian=TRUE]->cohortCounts
    
}

if("tabixResult" %in% names(opt)) {
    fread(opt$tabixResult)->exomeVariantTabixRes
    exomeVariantTabixRes[,AlleleKey:=paste(Chrom,Pos,Alternate,sep="_"),]
    exomeVariantTabixRes[,list(AlleleKey,Reference_EVS=Reference,MAF_EuropeanAmerican),]->exomeVariant_short
    setkey(exomeVariant_short,AlleleKey)
    
    exomeVariant_short[cohortCounts,allow.cartesian=TRUE]->cohortCounts
}

if("bamUnseqResult" %in% names(opt)) {
    fread(opt$bamUnseqResult)->bamUnseqs
    setnames(bamUnseqs,c('Chrom','Position'))
    setkeyv(bamUnseqs,c('Chrom','Position'))
    unique(bamUnseqs[,notSequencedIn:=.N,by=c("Chrom","Position")])->bamUnseqs
    
    setkeyv(cohortCounts,c('Chromosome','Position'))
    bamUnseqs[cohortCounts]->cohortCounts
    cohortCounts[is.na(notSequencedIn),notSequencedIn:=0,]
    cohortCounts[,Assumed_Total_Alleles:=Assumed_Total_Alleles-notSequencedIn,]
    cohortCounts[,Other_Alleles_Count:=Other_Alleles_Count-notSequencedIn,]
}

if("ref1000GPath" %in% names(opt)) {
  cohortCounts[,c("p.value_Xsqr_1000G","stat_Xsqr_1000G")
               :=runChiSqrd(Allele_Count,Other_Alleles_Count,
                            Allele_Count_1000G,Total_Alleles_1000G-Allele_Count_1000G
               ),by=AlleleKey]
  cohortCounts$adj_pvalue_bonferroni_1000G<-p.adjust(cohortCounts$p.value_Xsqr_1000G,method = "bonferroni")
  cohortCounts$adj_pvalue_fdr_1000G<-p.adjust(cohortCounts$p.value_Xsqr_1000G,method = "fdr")
  
}

# Chromosome    Position    Change    Cohort_Allele_Count    Cohort_Allele_Frequency    1000G_Allele_Count
# 1000G_Allele_Frequency    Chi-Squared_Statistic    P-Value    Bonferroni_Correction_0.05_Significance
# Gene    Effect    cDNA_pos    codon_pos    Protein_pos    AA_change    Grantham_score    dbsnp    PolyPhen
# SIFT    Protein_Domain    Clinic_Sig    Canonical_Trans    GERP    PHYLOP100
if("ref1000GPath" %in% names(opt) && "tabixResult" %in% names(opt)) {
  write.table(file = paste(opt$output,"significant1000G","withEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
   x = cohortCounts[adj_pvalue_fdr_1000G<0.05,list(
    Chrom,
    Position,
    Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
    Cohort_Allele_Count=Allele_Count,
    Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles,
    Allele_Count_1000G,
    Allele_Frequency_1000G=Allele_Count_1000G/Total_Alleles_1000G,
    stat_Xsqr_1000G,
    p.value_Xsqr_1000G,
    adj_pvalue_bonferroni_1000G,
    adj_pvalue_fdr_1000G,
    Gene, Effect, cDNA_pos, Codon_pos, Protein_pos, AA_change, 
    Grantham_Score, dbsnp, PolyPhen, SIFT, Protein_Domain, Clinic_Sig, Canonical_Trans, 
    GERP, PHYLOP100, MAF_EuropeanAmerican 
    ),][order(adj_pvalue_fdr_1000G)])
  write.table(file = paste(opt$output,"NotIn1000G","withEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[is.na(p.value_Xsqr_1000G),list(
                Chrom,
                Position,
                Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles,
                MAF_EuropeanAmerican 
              ),])
} else if(!("ref1000GPath" %in% names(opt)) && "tabixResult" %in% names(opt)) {
  write.table(file = paste(opt$output,"withEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[,list(
                Chrom,
                Position,
                Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles,
                MAF_EuropeanAmerican 
              ),])
} else if("ref1000GPath" %in% names(opt) && !("tabixResult" %in% names(opt))) {
  write.table(file = paste(opt$output,"significant1000G.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[adj_pvalue_fdr_1000G<0.05,list(
                Chrom,
                Position,
                Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles,
                Allele_Count_1000G,
                Allele_Frequency_1000G=Allele_Count_1000G/Total_Alleles_1000G,
                stat_Xsqr_1000G,
                p.value_Xsqr_1000G,
                adj_pvalue_bonferroni_1000G,
                adj_pvalue_fdr_1000G,
                Gene, Effect, cDNA_pos, Codon_pos, Protein_pos, AA_change, 
                Grantham_Score, dbsnp, PolyPhen, SIFT, Protein_Domain, Clinic_Sig, Canonical_Trans, 
                GERP, PHYLOP100
              ),][order(adj_pvalue_fdr_1000)])
  write.table(file = paste(opt$output,"NotIn1000G.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[is.na(p.value_Xsqr_1000G),list(
                Chrom,
                Position,
                Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles
              ),])
} else {
  write.table(file = paste(opt$output,"alone.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[,list(
                Chrom,
                Position,
                Change=paste("?->",Observed_Allele),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency=Allele_Count/Assumed_Total_Alleles,
              ),])
}
