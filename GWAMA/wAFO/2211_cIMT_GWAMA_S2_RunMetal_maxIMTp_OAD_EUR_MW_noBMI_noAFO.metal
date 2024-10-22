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


PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/QCed/OAD/maxRCF/EUR/CLEANED.CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922_3.cpaid_noDup_noAFO.tbl.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/NFBC/QCed/OAD/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.OAD.MW.noBMI.NFBC1966.EUR.JR.230801.cpaid_noDup_noAFO.tbl.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GRACE/QCed/OAD/maxLCF/EUR/CLEANED.GRACE_cIMT.maxLCF.OAD.MW.noBMI.EUR.AL.20240903_tmp.cpaid_noDup_noAFO.tbl.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/UKB/QCed/OAD/maxIMT240/BRIT/CLEANED.EGG_cIMT.OAD.maxIMT240_2.MW.noBMI.UKB.BRIT.AL.240920_tmp.cpaid_noDup_noAFO.tbl.gz  
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ALSPAC/QCed/OAD/maxLCF/EUR/CLEANED.ALSPAC_cIMT.maxLCF.OAD.MW.noBMI.EUR.AL.20240702.cpaid_noDup_noAFO.tbl.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/UKB/QCed/OAD/maxIMT240/NbEUR/CLEANED.EGG_cIMT.OAD.maxIMT240_2.MW.noBMI.UKB.NbEUR.AL.240920_tmp.cpaid_noDup_noAFO.tbl.gz
PROCESS /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/FAMILY/QCed/OAD/maxLCF/EUR/CLEANED.FAMILY_cIMT.maxLCF.OAD.MW.noBMI.EUR.AL.20240825_2.cpaid_noDup_noAFO.tbl.gz




OUTFILE /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAMA/OAD/maxLCFp/OAD_EUR/MW/noBMI/241015_maxLCFp_OAD_EUR_MW_noBMI_GWAMA_noAFO .TBL
ANALYZE HETEROGENEITY

