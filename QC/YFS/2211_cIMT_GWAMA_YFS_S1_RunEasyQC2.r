
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)
  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"

  qc2_p=paste0(m_p, "/Scripts/QC/YFS")
  o_p=paste0(m_p, "/Results/GWAS/YFS/QCed")

    # YFS YAD maxLCF  
    setwd(paste0(o_p, "/YAD/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/YFS_YAD_maxLCF_EUR_quant.ecf"))
     
    # YFS YAD maxLCF  
    setwd(paste0(o_p, "/ADL/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/YFS_ADL_maxLCF_EUR_quant.ecf"))
     
  
# d<-fread("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/YFS/Raw/EGG_cIMT.maxLCF.YAD.MW.adjBMI.YFS.EUR.LPL.231001.txt.gz")
