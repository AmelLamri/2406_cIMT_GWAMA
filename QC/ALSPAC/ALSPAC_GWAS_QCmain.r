
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf


    library(EasyQC2)
    library(data.table)

    m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
    qc2_p=paste0(m_p, "/Scripts/QC/ALSPAC/")
    o_p=paste0(m_p, "/Results/GWAS/ALSPAC/Cleaned")


    # YAD maxLCF
    setwd(paste0(o_p, "/YAD/maxLCF"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_YAD.ecf"))


    # OAD maxLCF
    setwd(paste0(o_p, "/OAD/maxLCF"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_OAD.ecf"))


    # ADL maxLCF
    setwd(paste0(o_p, "/ADL/maxLCF"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_ADL.ecf"))

    # AAD maxLCF
    setwd(paste0(o_p, "/AAD/maxLCF"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_AAD.ecf"))





