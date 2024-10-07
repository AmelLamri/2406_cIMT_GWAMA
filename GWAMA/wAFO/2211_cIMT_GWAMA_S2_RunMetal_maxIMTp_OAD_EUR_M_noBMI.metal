# Input columns:
MARKER cpaid
ALLELE EFFECT_ALLELE NON_EFFECT_ALLELE
EFFECT BETA
STDERRLABEL SE
FREQLABEL EAF
PVALUE PVAL
CUSTOMVARIABLE N
LABEL N AS N
SCHEME STDERR
WEIGHT N
#USESTRAND OFF
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
GENOMICCONTROL ON
#STRAND STRAND



PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/QCed/OAD/maxRCF/EUR/CLEANED.CP_cIMT.maxRCF.OAD.M.noBMI.EUR.KL.230922_2.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/NFBC/QCed/OAD/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.OAD.M.noBMI.NFBC1966.EUR.JR.230801.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GRACE/QCed/OAD/maxLCF/EUR/CLEANED.GRACE_cIMT.maxLCF.OAD.M.noBMI.EUR.AL.20240903_tmp.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/UKB/QCed/OAD/maxIMT240_2/BRIT/CLEANED.EGG_cIMT.OAD.maxIMT240_2.M.noBMI.UKB.BRIT.AL.240927_tmp.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ALSPAC/QCed/OAD/maxLCF/EUR/CLEANED.ALSPAC_cIMT.maxLCF.OAD.M.noBMI.EUR.AL.20240702.cpaid.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/UKB/QCed/OAD/maxIMT240_2/NbEUR/CLEANED.EGG_cIMT.OAD.maxIMT240_2.M.noBMI.UKB.NbEUR.AL.240920_tmp.cpaid.gz

OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/OAD/EUR/maxLCFp/noBMI/241005_maxIMT_OAD_EUR_M_noBMI_GWAMA .TBL
ANALYZE HETEROGENEITY

