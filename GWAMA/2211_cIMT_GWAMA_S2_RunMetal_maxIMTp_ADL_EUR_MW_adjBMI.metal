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


PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/QCed/ADL/maxRCF/EUR/CLEANED.CP_cIMT.maxRCF.ADL.MW.adjBMI.EUR.KL.230922_2.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/YFS/QCed/ADL/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.ADL.MW.adjBMI.YFS.EUR.LPL.231001.cpaid.gz     
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ALSPAC/QCed/ADL/maxLCF/EUR/CLEANED.ALSPAC_cIMT.maxLCF.ADL.MW.adjBMI.EUR.AL.20240702.cpaid.gz

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/ADL/EUR/maxLCFp/adjBMI/240926_maxIMT_ADL_EUR_MW_adjBMI_GWAMA_ .TBL
ANALYZE HETEROGENEITY

