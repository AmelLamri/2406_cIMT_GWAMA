
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf

DEFINE	--pathOut /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/GRACE/Cleaned/EUR/maxLCF/


        library(EasyQC2)
        library(data.table)
        m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
        qc2_p=paste0(m_p, "/Scripts/QC/")
        o_p=paste0(m_p, "/Results/GWAS")

    # GRACE EUR 
        setwd(paste0(o_p, "/GRACE/Cleaned/EUR/maxLCF"))
    EasyQC2(paste0(qc2_p, "/GRACE/EGG_cIMT_EasyQC2_GRACE_EUR.ecf"))


    # GRACE AMR 
        setwd(paste0(o_p, "/GRACE/Cleaned/AMR/maxLCF"))
    EasyQC2(paste0(qc2_p, "/GRACE/EGG_cIMT_EasyQC2_GRACE_AMR.ecf"))

