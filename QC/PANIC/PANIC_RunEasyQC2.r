
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/PANIC")
  ir1_p=paste0(m_p, "/Results/GWAS/PANIC/")
  o1_p=paste0(ir1_p, "/QCed/")
 

 # Run EasyQC 
  # ADO maxLCF EUR 
    setwd(paste0(o1_p, "/ADO/meanLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/PANIC_ADO_meanLCF_EUR_QC_quant.ecf"))
