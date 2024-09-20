
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf


  library(EasyQC2)
  library(data.table)
  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"

  qc2_p=paste0(m_p, "/Scripts/QC/FAMILY")
  o_p=paste0(m_p, "/Results/GWAS")

  setwd(paste0(o_p, "/FAMILY/Cleaned/ADO/maxLCF/"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_ADO_FAMILY_QC_quant.ecf"))

  setwd(paste0(o_p, "/FAMILY/Cleaned/ADL/maxLCF/"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_ADL_FAMILY_QC_quant.ecf"))

  setwd(paste0(o_p, "/FAMILY/Cleaned/OAD/maxLCF/"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_OAD_FAMILY_QC_quant.ecf"))

  setwd(paste0(o_p, "/FAMILY/Cleaned/AAD/maxLCF/"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_AAD_FAMILY_QC_quant.ecf"))





