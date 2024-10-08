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


PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ABCD/QCed/ADO/maxCCA/EUR/CLEANED.EGG_cIMT.maxCCA.ADO.MW.adjBMI.EUR.EM.20230626.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/QCed/ADO/maxRCF/EUR/CLEANED.CP_cIMT.maxRCF.ADO.MW.adjBMI.EUR.KL.230922_2.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/FAMILY/QCed/ADO/maxLCF/EUR/CLEANED.FAMILY_cIMT.maxLCF.ADO.MW.adjBMI.EUR.AL.20240825.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR1/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.230731.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR2/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.231005.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/PANIC/QCed/ADO/meanLCF/EUR/CLEANED.PANIC_cIMT.meanLCF.ADO.MW.adjBMI.EUR.TS.231213.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/SWS/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.adjBMI.EUR.PT.240130.cpaid.gz



ANALYZE HETEROGENEITY

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/ADO/EUR/maxLCFp/adjBMI/240926_maxIMT_ADO_EUR_MW_adjBMI_GWAMA_ .TBL
