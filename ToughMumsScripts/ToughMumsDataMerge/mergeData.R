suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

options_list<-list(
  make_option(c("--cohortCounts","-c"),help="cohorts count file"),
  make_option(c("--output","-o"),help="ouput full path name"),
  make_option(c("--tabixResult","-t"),help="full path to tabix processed result"),
  make_option(c("--femaleCount","-f"),help="Number of females",type="integer"),
  make_option(c("--maleCount","-m"),help="Number of males",type="integer"),
  make_option(c("--ref1000GPath","-r"),help="Path to 1000 Genomes reference file"),
  make_option(c("--bamUnseqResult","-u"),help="Path to file with the file names containing the positions that are not sequenced per sample")
)

opt<-parse_args(OptionParser(option_list = options_list),positional_arguments = FALSE)

t0<-Sys.time()

# we index all tables by a new column: chr-position-allele
# cohort columns: Chromosome      Position        Observed_Allele Allele_Count    Assumed_Total_Alleles
fread(opt$cohortCounts)->cohortCounts
cohortCounts[,Assumed_Total_Alleles:=2*(opt$femaleCount+opt$maleCount),]
cohortCounts[Chromosome=="Y",Assumed_Total_Alleles:=opt$maleCount,]
cohortCounts[,Other_Alleles_Count:=Assumed_Total_Alleles-Allele_Count,]
cohortCounts[,AlleleKey:=paste(Chromosome,Position,Observed_Allele,sep="_"),]
setkey(cohortCounts,AlleleKey)

t1<-Sys.time()
print(paste("Initial proc of cohortCounts :",t1-t0,sep=""))

runChiSqrd<-function(cohortAllele,cohortOther,referenceAllele,referenceOther) {
  if(is.na(referenceAllele) || referenceAllele < 0 
     || is.na(referenceOther) || referenceOther < 0
     || is.na(cohortAllele) || cohortAllele < 0 
     || is.na(cohortOther)  || cohortOther < 0) {
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
    setnames(reference,
             old=c('cDNA pos','codon pos','Protein pos','AA change','Grantham score','Protein Domain','Clinic Sig','Canonical Trans'),
             new=c('cDNA_pos','Codon_pos','Protein_pos','AA_change','Grantham_score','Protein_Domain','Clinic_Sig','Canonical_Trans'))
    
    reference[,list(
      Obs_Allele_1000G=Observed_Allele[1],
      Allele_Count_1000G=Allele_Count[1],
      Total_Alleles_1000G=Total_Alleles[1],
      Gene=Gene[1],Ref_Allele_1000G=Ref_Allele[1],
      cDNA_pos=cDNA_pos[1],
      Codon_pos=Codon_pos[1],
      Protein_pos=Protein_pos[1],
      AA_change=AA_change[1],
      Grantham_Score=Grantham_score[1],
      dbsnp=paste(unique(unlist(strsplit(dbsnp,split = ";"))),collapse = ";"), 
      PolyPhen=PolyPhen[1],SIFT=SIFT[1],
      Protein_Domain=paste(unique(unlist(strsplit(Protein_Domain,split='&'))),collapse = "&"),
      Clinic_Sig=paste(unique(Clinic_Sig),collapse="-"),
      Canonical_Trans=gsub(pattern = '(^-|-$)', replacement = "", x = paste(unique(Canonical_Trans),collapse="-")),
      GERP=GERP[1],PHYLOP100=PHYLOP100[1]),
    by=.(AlleleKey,Effect)]->ref_short
    
    setkey(ref_short,AlleleKey)
    
    ref_short[cohortCounts,allow.cartesian=TRUE]->cohortCounts
    
}

t2<-Sys.time()
print(paste("Initial processing of 1000G : ",t2-t1,sep=""))

if("tabixResult" %in% names(opt)) {
    fread(opt$tabixResult)->exomeVariantTabixRes
    exomeVariantTabixRes[,EA_Alt_Count:=as.numeric(EA_Alt_Count),]
    exomeVariantTabixRes[,MAF_EuropeanAmerican:=as.numeric(MAF_EuropeanAmerican),]
    exomeVariantTabixRes[,AlleleKey:=paste(Chrom,Pos,Alternate,sep="_"),]
    exomeVariantTabixRes[MAF_EuropeanAmerican>0,Allele_Count_EA_Others_EVS:=EA_Alt_Count/MAF_EuropeanAmerican,]
    exomeVariantTabixRes[,list(AlleleKey,Reference_EVS=Reference,MAF_EuropeanAmerican,Allele_Count_EA_EVS=EA_Alt_Count,Allele_Count_EA_Others_EVS),]->exomeVariant_short
    setkey(exomeVariant_short,AlleleKey)
    
    exomeVariant_short[cohortCounts,allow.cartesian=TRUE]->cohortCounts
}

t3<-Sys.time()
print(paste("Initial processing of EVS : ",t3-t2,sep=""))

if("bamUnseqResult" %in% names(opt)) {
    fread(opt$bamUnseqResult)->bamUnseqs
    setnames(bamUnseqs,c('Chrom','Position'))
    setkeyv(bamUnseqs,c('Chrom','Position'))
    unique(bamUnseqs[,notSequencedIn:=.N,by=c("Chrom","Position")])->bamUnseqs
    
    setkeyv(cohortCounts,c('Chromosome','Position'))
    bamUnseqs[cohortCounts]->cohortCounts
    cohortCounts[is.na(notSequencedIn),notSequencedIn:=0,]
    cohortCounts[,Assumed_Total_Alleles:=Assumed_Total_Alleles-2*notSequencedIn,]
    cohortCounts[,Other_Alleles_Count:=Other_Alleles_Count-2*notSequencedIn,]
}

t4<-Sys.time()
print(paste("Initial processing of BAMs unseq : ",t4-t3,sep=""))

cohortCounts[,Cohort_Allele_Frequency:=Allele_Count/Assumed_Total_Alleles,]

setkey(cohortCounts,AlleleKey)
if("ref1000GPath" %in% names(opt)) {
  cohortCounts[,c("p.value_Xsqr_1000G","stat_Xsqr_1000G")
               :=runChiSqrd(Allele_Count,Other_Alleles_Count,
                            Allele_Count_1000G,Total_Alleles_1000G-Allele_Count_1000G
               ),by=AlleleKey]
  unique(cohortCounts[,.(Allele_Count,Other_Alleles_Count,
                  Allele_Count_1000G,Total_Alleles_1000G),by=AlleleKey])->uniqueForPValues
  setkey(uniqueForPValues,AlleleKey)
  uniqueForPValues[,c("p.value_Xsqr_1000G","stat_Xsqr_1000G")
                   :=runChiSqrd(Allele_Count,Other_Alleles_Count,
                                Allele_Count_1000G,Total_Alleles_1000G-Allele_Count_1000G
                   ),by=AlleleKey]
  uniqueForPValues[,adj_pvalue_bonferroni_1000G:=p.adjust(p.value_Xsqr_1000G,method = "bonferroni"),]
  uniqueForPValues[,adj_pvalue_fdr_1000G:=p.adjust(p.value_Xsqr_1000G,method = "fdr"),]
  uniqueForPValues[,c('Allele_Count','Other_Alleles_Count','Allele_Count_1000G','Total_Alleles_1000G'):=NULL,]
  uniqueForPValues[cohortCounts]->cohortCounts
  cohortCounts[Total_Alleles_1000G>0,Allele_Frequency_1000G:=Allele_Count_1000G/Total_Alleles_1000G,]
  cohortCounts[!is.na(Allele_Frequency_1000G),OverAbundance_Cohort_1000G:=Cohort_Allele_Frequency/Allele_Frequency_1000G,]
}

t5<-Sys.time()
print(paste("Secondary proc of 1000G : ",t5-t4,sep=""))

if("tabixResult" %in% names(opt)) {
  
  unique(cohortCounts[,.(Allele_Count,Other_Alleles_Count,
                         Allele_Count_EA_EVS,Allele_Count_EA_Others_EVS),by=AlleleKey])->uniqueForPValues
  setkey(uniqueForPValues,AlleleKey)
  uniqueForPValues[,c("p.value_Xsqr_EA_EVS","stat_Xsqr_EA_EVS")
                   :=runChiSqrd(Allele_Count,Other_Alleles_Count,
                                Allele_Count_EA_EVS,Allele_Count_EA_Others_EVS
                   ),by=AlleleKey]
  uniqueForPValues[,adj_pvalue_bonferroni_EA_EVS:=p.adjust(p.value_Xsqr_EA_EVS,method = "bonferroni"),]
  uniqueForPValues[,adj_pvalue_fdr_EA_EVS:=p.adjust(p.value_Xsqr_EA_EVS,method = "fdr"),]
  uniqueForPValues[,c('Allele_Count','Other_Alleles_Count','Allele_Count_EA_EVS','Allele_Count_EA_Others_EVS'):=NULL,]
  uniqueForPValues[cohortCounts]->cohortCounts
  cohortCounts[!is.na(Allele_Count_EA_Others_EVS),Allele_Frequency_EA_EVS:=Allele_Count_EA_EVS/(Allele_Count_EA_EVS+Allele_Count_EA_Others_EVS),]
  cohortCounts[!is.na(Allele_Frequency_EA_EVS) && Allele_Frequency_EA_EVS>0,OverAbundance_Cohort_EA_EVS:=Cohort_Allele_Frequency/Allele_Frequency_EA_EVS,]
  
}

t6<-Sys.time()
print(paste("Secondary processing of EVS : ",t6-t5,sep=""))

# Chromosome    Position    Change    Cohort_Allele_Count    Cohort_Allele_Frequency    1000G_Allele_Count
# 1000G_Allele_Frequency    Chi-Squared_Statistic    P-Value    Bonferroni_Correction_0.05_Significance
# Gene    Effect    cDNA_pos    codon_pos    Protein_pos    AA_change    Grantham_score    dbsnp    PolyPhen
# SIFT    Protein_Domain    Clinic_Sig    Canonical_Trans    GERP    PHYLOP100
if("ref1000GPath" %in% names(opt) && "tabixResult" %in% names(opt)) {
  write.table(file = paste(opt$output,"significant1000G","withEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
   x = cohortCounts[adj_pvalue_fdr_1000G<0.05 | adj_pvalue_fdr_EA_EVS<0.05,list(
    Chrom,
    Position,
    Change=paste(Reference_Allele,Observed_Allele,sep="->"),
    Cohort_Allele_Count=Allele_Count,
    Cohort_Allele_Frequency,
    Allele_Count_1000G,
    Allele_Frequency_1000G,
    stat_Xsqr_1000G,
    p.value_Xsqr_1000G,
    adj_pvalue_bonferroni_1000G,
    adj_pvalue_fdr_1000G,
    OverAbundance_Cohort_1000G,
    Allele_Count_EA_EVS,
    Allele_Frequency_EA_EVS,
    stat_Xsqr_EA_EVS,
    p.value_Xsqr_EA_EVS,
    adj_pvalue_bonferroni_EA_EVS,
    adj_pvalue_fdr_EA_EVS,
    OverAbundance_Cohort_EA_EVS,
    Gene, Effect, cDNA_pos, Codon_pos, Protein_pos, AA_change, 
    Grantham_Score, dbsnp, PolyPhen, SIFT, Protein_Domain, Clinic_Sig, Canonical_Trans, 
    GERP, PHYLOP100, MAF_EuropeanAmerican 
    ),][order(adj_pvalue_fdr_1000G),])
  
  cohortCounts[is.na(p.value_Xsqr_1000G) & is.na(p.value_Xsqr_EA_EVS),list(
    Chrom,
    Position,
    Change=paste(Reference_Allele,Observed_Allele,sep="->"),
    Cohort_Allele_Count=Allele_Count,
    Cohort_Allele_Frequency#,
    #Log=log2(Cohort_Allele_Frequency)
  ),]->toWrite
  #toWrite[!is.infinite(Log),LogNormalPValue:=as.vector(2*pnorm(-abs(scale(Log,center = TRUE,scale = TRUE)))),]
  write.table(file = paste(opt$output,"NotIn1000G","NotInEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
             #x = toWrite[order(LogNormalPValue),])
             x = toWrite[order(-Cohort_Allele_Frequency),])
} else if(!("ref1000GPath" %in% names(opt)) && "tabixResult" %in% names(opt)) {
  write.table(file = paste(opt$output,"withEVS.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[,list(
                Chrom,
                Position,
                Change=paste(Reference_Allele,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency,
                Allele_Count_EA_EVS,
                Allele_Frequency_EA_EVS,
                stat_Xsqr_EA_EVS,
                p.value_Xsqr_EA_EVS,
                adj_pvalue_bonferroni_EA_EVS,
                adj_pvalue_fdr_EA_EVS,
                OverAbundance_Cohort_EA_EVS,
                MAF_EuropeanAmerican 
              ),])
} else if("ref1000GPath" %in% names(opt) && !("tabixResult" %in% names(opt))) {
  write.table(file = paste(opt$output,"significant1000G.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[adj_pvalue_fdr_1000G<0.05,list(
                Chrom,
                Position,
                Change=paste(Ref_Allele_1000G,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency,
                Allele_Count_1000G,
                Allele_Frequency_1000G,
                stat_Xsqr_1000G,
                p.value_Xsqr_1000G,
                adj_pvalue_bonferroni_1000G,
                adj_pvalue_fdr_1000G,
                OverAbundance_Cohort_1000G,
                Gene, Effect, cDNA_pos, Codon_pos, Protein_pos, AA_change, 
                Grantham_Score, dbsnp, PolyPhen, SIFT, Protein_Domain, Clinic_Sig, Canonical_Trans, 
                GERP, PHYLOP100
              ),][order(adj_pvalue_fdr_1000)])
  write.table(file = paste(opt$output,"NotIn1000G.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[is.na(p.value_Xsqr_1000G),list(
                Chrom,
                Position,
                Change=paste(Reference_Allele,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency
              ),])
} else {
  write.table(file = paste(opt$output,"alone.xls",sep="_"),sep = '\t', row.names = F, quote = F, 
              x = cohortCounts[,list(
                Chrom,
                Position,
                Change=paste(Reference_Allele,Observed_Allele,sep="->"),
                Cohort_Allele_Count=Allele_Count,
                Cohort_Allele_Frequency
              ),])
}
