# Input columns:
MARKER cpaid
ALLELE EFFECT_ALLELE NON_EFFECT_ALLELE
EFFECT BETA
STDERRLABEL SE
FREQLABEL EAF
PVALUE PVAL
STRAND STRAND
CUSTOMVARIABLE N
LABEL N AS N
SCHEME STDERR
WEIGHT N
USESTRAND ON
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
GENOMICCONTROL ON

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch2/Cleaned/CLD_LCF/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.231005.cpaid.1_clean.gz 

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch1/Cleaned/CLD_LCF/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.230912.cpaid_clean.gz

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Cleaned/ADO_CCA/CLEANED.EGG_cIMT.meanCCA.ADO.MW.noBMI.EUR.EM.20230626.cpaid.gz

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Cleaned/ADO/maxRCF/CLEANED.CP_cIMT.maxRCF.ADO.MW.noBMI.EUR.KL.230922.cpaid.gz


PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/PANIC/Cleaned/ADO/meanLCF/CLEANED.PANIC_cIMT.meanLCF.ADO.MW.adjBMI.EUR.TS.231213.cpaid.gz

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/FAMILY/Cleaned/ADO/maxLCF/CLEANED.EGG_cIMT.maxIMT.ADO.MW.noBMI.EUR.AL.20240528.cpaid.gz

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/ADO/maxIMT/EGG_ADO_MaxIMT_MW_noBMI_AL_240909.TBL


ANALYZE HETEROGENEITY

# TO DO add CHCP non Euro  
