
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/ABCD")
  ir_p=paste0(m_p, "/Results/GWAS/ABCD/")
  rw_p=paste0(ir_p, "/Raw/")
  md_p=paste0(ir_p, "/Modified/")
  o_p=paste0(ir_p, "/QCed")
# Run EasyQC 
  # ADO maxCCA EUR  
    setwd(paste0(o_p, "/ADO/maxCCA/EUR"))
    EasyQC2(paste0(qc2_p, "/ABCD_ADO_maxCCA_EUR_QC_quant.ecf"))
