
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



        library(EasyQC2)
        library(data.table)
        m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
        qc2_p=paste0(m_p, "/Scripts/QC/")
        o_p=paste0(m_p, "/Results/GWAS")

    # UKB BRIT 
    setwd(paste0(o_p, "/UKB/Cleaned/BRIT/maxIMT240_2/"))
    EasyQC2(paste0(qc2_p, "/UKB/EGG_cIMT_EasyQC2_UKB_BRIT.ecf"))




