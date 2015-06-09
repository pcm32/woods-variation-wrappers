
# on this location the vcf files should be stored, either combined or separated
# SNPs and indels.
SHORTQUEUE=short
LONGQUEUE=long

WRAPPERDIR=/data/woods/variation-wrappers

OLD_VCFS_PATH=/data/woods/Old_VCFs
VEP_DIR=$OLD_VCFS_PATH/variant_effect_predictor
FORMATTING_TOOLS_PATH=$OLD_VCFS_PATH/Formatting_Tools
FORMATTINGTEMP=/results/woods/DataFormattingTemp
ANNOTATED_VCFS_PATH=/data/woods/Annotated_VCFs

GENEFILTERPATH=/data/woods/GeneFilter
GENEFILTERRESPATH=/results/woods/GeneFilter

MUTATIONFILTERPATH=/data/woods/MutationFilter
MUTATIONFILTERRESULTS=/results/woods/MutationFilter
MUTFILTXREFsDIR=$ANNOTATED_VCFS_PATH/Cross_Reference

RUNJOBSPATH=/data/woods/runJobs
CASECONTROLRESPATH=/results/woods/CaseControl
CASECONTROLPATH=/data/woods/CaseControl
CASECONTROLTEMP=/results/woods/CaseControl

BAMSPATH=/data/woods/Bam_Files
TOUGHMUMSRESULTS=/results/woods/ToughMums
TOUGHMUMSPATH=/data/woods/ToughMums
TOUGHMUMSPATHSC=$WRAPPERDIR/ToughMumsScripts

PATH=/shared/software/bedtools-2.24.0/bin:$PATH
