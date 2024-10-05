
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/YAD/")
  o_p=paste0(m_p, "/Results/GWAS")

    # YFS 
    setwd(paste0(o_p, "//YFS/Cleaned/AAD_maxLCF"))
    EasyQC2(paste0(qc2_p, "/2211_cIMT_YAD_GWAMA_S1_RunEasyQC2_1_YFS.ecf"))
     
    # ALSPAC
    setwd(paste0(o_p, "//GenR/ALSPAC/Cleaned/AAD_maxLCF"))
    EasyQC2(paste0(qc2_p, "/"))
