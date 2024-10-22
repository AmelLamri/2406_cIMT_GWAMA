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

PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/YFS/QCed/YAD/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.YAD.MW.adjBMI.YFS.EUR.LPL.231001.cpaid_noDup_noAFO.tbl.gz  
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ALSPAC/QCed/YAD/maxLCF/EUR/CLEANED.ALSPAC_cIMT.maxLCF.YAD.MW.adjBMI.EUR.AL.20240702.cpaid_noDup_noAFO.tbl.gz

ANALYZE HETEROGENEITY

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAMA/YAD/maxLCFp/YAD_EUR/MW/adjBMI/241015_maxLCFp_YAD_EUR_MW_adjBMI_GWAMA_noAFO .TBL
