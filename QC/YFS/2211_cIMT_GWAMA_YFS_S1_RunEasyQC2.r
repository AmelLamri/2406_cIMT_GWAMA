
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)
  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"

  qc2_p=paste0(m_p, "/Scripts/QC/")
  o_p=paste0(m_p, "/Results/GWAS")

    # YFS YAD maxLCF  
    setwd(paste0(o_p, "/YFS/Cleaned/YAD/maxLCF"))
    EasyQC2(paste0(qc2_p, "/YFS/2211_cIMT_YFS_YAD_maxLCF_S2_RunEasy.ecf"))
     
    # YFS YAD maxLCF  
    setwd(paste0(o_p, "/YFS/Cleaned/ADL/maxLCF"))
    EasyQC2(paste0(qc2_p, "/YFS/2211_cIMT_YFS_ADL_maxLCF_S2_RunEasy.ecf"))
     