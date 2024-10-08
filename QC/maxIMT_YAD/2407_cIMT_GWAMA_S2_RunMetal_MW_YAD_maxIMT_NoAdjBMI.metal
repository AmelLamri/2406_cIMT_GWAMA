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

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/YFS/Cleaned/YAD_maxLCF/CLEANED.EGG_cIMT.maxLCF.YAD.MW.noBMI.YFS.EUR.LPL.231001.cpaid.gz 
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ALSPAC/Cleaned/YAD_maxLCF/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.230912.cpaid_clean.gz

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/YAD_maxLCF/EGG_YAD_MaxIMT_MW_noBMI_AL_240704.TBL





ANALYZE HETEROGENEITY



