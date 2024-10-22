
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



        library(EasyQC2)
        library(data.table)
        m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
        qc2_p=paste0(m_p, "/Scripts/QC/")
        o_p=paste0(m_p, "/Results/GWAS")

    # OAD BRIT 
    setwd(paste0(o_p, "/UKB/QCed/OAD/maxIMT240_2/BRIT"))
    EasyQC2(paste0(qc2_p, "/UKB/EGG_cIMT_EasyQC2_UKB_OAD_BRIT_max240_2.ecf"))

    
    # OAD NBEUR 
    setwd(paste0(o_p, "/UKB/QCed/OAD/maxIMT240_2/NBEUR"))
    EasyQC2(paste0(qc2_p, "/UKB/EGG_cIMT_EasyQC2_UKB_OAD_NBEUR_max240_2.ecf"))


    