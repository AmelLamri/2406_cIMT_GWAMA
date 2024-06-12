
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)
    library(tidyr)

    # GenR 
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2/GenR/Batch1/EGG_cIMT_GenR_B1_QC_quant.ecf")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2/GenR/Batch2/EGG_cIMT_GenR_B2_QC_quant.ecf")

    # YFS 
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_YFS_QC_quant.ecf")
     
    # NFBC1966
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_NFBC1966_QC_quant.ecf")

    # CHCP 
	EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_CHCP_QC_quant.ecf")
    # embedded nul(s) found in input
    # embedded nul(s) found in input
    # EASY WARNING:
    #  There is no TAB in the header of file
    # /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw_GWAS/CP_cIMT.maxRCFxAGE.MW.adjBMI.EUR.KL.230922.zip
    #  Please make sure that you have defined the correct delimiter!
    # Error in GWADATA.init(container.GWADATA[[i]]) : EASY ERROR:GWADATA
    #  Defined column
    # CHR;POS;INFO;EFFECT_ALLELE;NON_EFFECT_ALLELE;EAF;N;BETA;SE;PVAL not available in file
    # /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw_GWAS/CP_cIMT.maxRCFxAGE.MW.adjBMI.EUR.KL.230922.zip
    #  Please check !!!
    # 
    # [1] FALSE
    # Warning messages:
    # 1: In scan(file = object@fileIn, what = "character", skip = numRowSkip,  :     
    #   embedded nul(s) found in input
    # 2: In scan(file = object@fileIn, what = "character", skip = numRowSkip,  :     
    #   embedded nul(s) found in input
    # 3: In GWADATA.init(container.GWADATA[[i]]) : EASY WARNING:
    #  There is no TAB in the header of file
    # /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw_GWAS/CP_cIMT.maxRCFxAGE.MW.adjBMI.EUR.KL.230922.zip
    #  Please make sure that you have defined the correct delimiter!
	
    # PANIC 
	EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_PANIC_QC_quant.ecf")
    # Error in GWADATA.init(container.GWADATA[[i]]) : EASY ERROR:GWADATA
    # Defined column CHR;POS;INFO;EFFECT_ALLELE;NON_EFFECT_ALLELE;EAF;N;BETA;SE;PVAL not available in file
    # /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/PANIC/Raw/PANIC_cIMT.meanLCF.ADO.M.adjBMI.EUR.TS.231213.txt.gz
    # Please check !!!
	
    # ABCD
    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Cleaned/")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_ABCD_QC_quant.ecf")
    #   +  EGG_cIMT.maxCCA.ADO.MEN.adjBMI.EUR.EM.20230626.txt.gz ->
    # Error in value[[3L]](cond) : line 1 did not have 19 elements
    # EASY ERROR:
    # Cannot read 'line 1 did not have 19 elements' from row 'NA' !!!
    # Please specifiy correct column class in --acolInClasses .

    # PANIC 
    
    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/PANIC/Cleaned/")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_PANIC_QC_quant.ecf")

    # CHCP  

    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Cleaned/maxRCF_ADO/EUR")
    
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//CHCP/EGG_cIMT_CHCP_maxRCF_ADO_EUR_QC_quant.ecf")

    
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_CHCP_meanRCF_QC_quant.ecf")



    # GenR Batch1 

    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch1/Cleaned/")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_maxLCF_GenR_B1_QC_quant.ecf")

    # GenR Batch 2 

    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch2/Cleaned/")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_maxLCF_GenR_B2_QC_quant.ecf")



 # Family  

    setwd("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/FAMILY/Cleaned/")
    EasyQC2("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts/EasyQC2//EGG_cIMT_maxLCF_ADO_FAMILY_QC_quant.ecf")




sudo mv /home/lamria/nutri/FAMILY/Projects/2301_cIMT_GWAS/Outputs/*.assoc.linear /hightide/lamria/nutri/FAMILY/Projects/2301_cIMT_GWAS/Outputs


genb2<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch2/Raw/EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.231005.gz"), nrow=100, head=T)
genb1<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Batch1/Raw/EGG_cIMT.maxLCF.CLD.MEN.adjBMI.GenR.EUR.MMG.230912.gz"), nrow=100, head=T)
yfs<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GenR/Raw/EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.231005.gz"), nrow=100, head=T)
nfbc<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/NFBC1966/Raw/EGG_cIMT.maxLCF.OAD.MEN.adjBMI.NFqBC1966.EUR.JR.230807.txt.gz"), nrow=100, head=T)
nfbc<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/NFBC1966/Raw/EGG_cIMT.maxLCF.OAD.MEN.adjBMI.NFqBC1966.EUR.JR.230807.txt.gz"), nrow=100, head=T)
fam<- read.table(gzfile ("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/FAMILY/Raw/EGG_cIMT.maxIMT.ADO.MW.noBMI.EUR.AL.20240528.txt.gz"), nrow=100, head=T)
chcp <-  read.table(unz("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw/CP_cIMT.maxRCF.ADO.MEN.adjBMI.EUR.KL.230922.zip"), nrow=100, head=T)
 colnames(nfbc)
#  [1] "SNPID"             "STRAND"            "BUILD"
#  [4] "CHR"               "POS"               "EFFECT_ALLELE"     
#  [7] "NON_EFFECT_ALLELE" "N"                 "N0"
# [10] "N1"                "N2"                "HWE_P"
# [13] "CALL_RATE"         "BETA"              "SE"
# [16] "PVAL"              "INFO_TYPE"         "INFO"



nfbc$EAF <- (nfbc$N0 + nfbc$N1)/  (nfbc$N0 + nfbc$N1+ nfbc$N2)

#CHCP 

#MEN have a lot of NA in EAF, so check why 

chcp <-  read.delim(unz("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw/CP_cIMT.maxRCF.ADO.MEN.adjBMI.EUR.KL.230922.zip", "CP_cIMT.maxRCF.ADO.MEN.adjBMI.EUR.KL.230922.txt"), head=T)
