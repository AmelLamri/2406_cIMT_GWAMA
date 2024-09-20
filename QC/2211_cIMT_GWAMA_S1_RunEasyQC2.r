
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)
  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"

  qc2_p=paste0(m_p, "/Scripts/EasyQC2/")
  o_p=paste0(m_p, "/Results/GWAS")

    # YFS 
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_YFS_QC_quant.ecf"))
     
    # NFBC1966
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_NFBC1966_QC_quant.ecf"))

    # CHCP 
	EasyQC2(paste0(qc2_p, "/EGG_cIMT_CHCP_QC_quant.ecf"))
   	
   
    # ABCD
    setwd(paste0(o_p, "/ABCD/Cleaned/ADO_CCA"))
    EasyQC2(paste0(qc2_p, "/ABCD/ADO_CCA/EGG_cIMT_ABCD_ADO_meanCCA_QC_quant.ecf"))
    EasyQC2(paste0(qc2_p, "/ABCD/ADO_CCA/EGG_cIMT_ABCD_ADO_maxCCA_QC_quant.ecf"))


    # PANIC 
    panicp="ADO_meanLCF" 
    setwd(paste0(o_p, "//PANIC/Cleaned/ADO_meanLCF"))
    EasyQC2(paste0(qc2_p, "/PANIC/EGG_cIMT_PANIC_ADO_meanLCF_QC_quant.ecf"))

    # CHCP  
    setwd(paste0(o_p, "//CHCP/Cleaned/maxRCF_ADO/EUR")
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_ADO_EUR_QC_quant.ecf"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_CHCP_meanRCF_QC_quant.ecf"))



    # GenR Batch1 

    setwd(paste0(o_p, "//GenR/Batch1/Cleaned/CLD_LCF"))
    EasyQC2(paste0(qc2_p, "/GenR/Batch1/EGG_cIMT_CLD_maxLCF_GenR_B1_QC_quant.ecf"))

    # GenR Batch 2 

    setwd(paste0(o_p, "//GenR/Batch2/Cleaned/CLD_LCF"))
    EasyQC2(paste0(qc2_p, "/GenR/Batch2/EGG_cIMT_CLD_maxLCF_GenR_B2_QC_quant.ecf"))



 # Family  

    setwd(paste0(o_p, "//FAMILY/Cleaned/"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_ADO_FAMILY_QC_quant.ecf"))




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


explore output 

    genR <- read.table(gzfile(paste0(o_p, "//GenR/Batch1/Cleaned/CLD_LCF/CLEANED.EGG_cIMT.maxLCF.CLD.MEN.adjBMI.GenR.EUR.MMG.230912.cpaid.gz")), head=T, nrow=1000)


#Get Ns 
  # YFS 
  
    #YAD 
    
      yfs_yad <- read.table(gzfile(paste0(o_p, "/YFS/Raw/EGG_cIMT.maxLCF.YAD.WOMEN.adjBMI.YFS.EUR.LPL.231001.txt.gz")), nrow=100, head=T)
    
    
    # ADL 
      yfs_adl <- read.table(gzfile(paste0(o_p, "/YFS/Raw/EGG_cIMT.maxLCF.ADL.MW.adjBMI.YFS.EUR.LPL.231001.txt.gz")), nrow=100, head=T)

  # FAM 
    #ADL 
      fam_adl <- read.table(gzfile(paste0(o_p, "/FAM/")), nrow=100, head=T)

EGG_cIMT.maxLCF.ADL.MW.adjBMI.YFS.EUR.LPL.231001.txt.gz


  # CHCP 
    # ADL 

      chcp_adl <- read.table(unz(paste0(o_p, "/YFS/Raw/EGG_cIMT.maxLCF.YAD.WOMEN.adjBMI.YFS.EUR.LPL.231001.txt.gz")), nrow=100, head=T)
      cat "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/Raw/CP_cIMT.minRCN.OAD.MW.noBMI.EUR.KL.230922 - Copy/CP_cIMT.minRCN.OAD.MW.noBMI.EUR.KL.230922.txt" | head 

    #OAD 

      chcp_adl <- read.table(paste0(o_p, "/CHCP/Raw/CP_cIMT.minRCN.OAD.MW.noBMI.EUR.KL.230922.zip"), nrow=100, skip=10)
      cat "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/Raw/CP_cIMT.minRCN.OAD.MW.noBMI.EUR.KL.230922 - Copy/CP_cIMT.minRCN.OAD.MW.noBMI.EUR.KL.230922.txt" | head 

      cat "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/CHCP/Raw/CP_cIMT.maxRCF.AAD.MW.adjBMI.EUR.KL.230922/CP_cIMT.maxRCF.AAD.MW.adjBMI.EUR.KL.230922.txt" | head 

      CP_cIMT.maxRCF.AAD.MW.adjBMI.EUR.KL.230922


   #NFBC 
     # OAD 
       nfbc_oad <- read.table(gzfile(paste0(o_p, "/NFBC1966/Raw/EGG_cIMT.maxLCF.OAD.MW.adjBMI.NFBC1966.EUR.JR.230807.txt.gz")), nrow=100, head=T)

  #UKB 

  cat "/genetics/Users/lamria/UKB/2405_cIMT_GWAS/Outputs/BRIT/maxIMT120_2/EGG_cIMT.OAD.maxIMT120_2.MEN..UKB.BRIT.AL.240502.1.out" | head -n 20
  

