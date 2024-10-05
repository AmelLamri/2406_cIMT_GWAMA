
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/GenR")
  ir1_p=paste0(m_p, "/Results/GWAS/GenR1/")
  ir2_p=paste0(m_p, "/Results/GWAS/GenR2/")
  o1_p=paste0(ir1_p, "/QCed/")
  o2_p=paste0(ir2_p, "/QCed/")

 # Run EasyQC 
  #  GenR1 ADO maxLCF EUR 
    setwd(paste0(o1_p, "/ADO/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/GenR1_ADO_maxLCF_EUR_QC_quant.ecf"))


    
  #  GenR2 ADO maxLCF EUR 
    setwd(paste0(o2_p, "/ADO/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/GenR2_ADO_maxLCF_EUR_QC_quant.ecf"))
    
# data has embeded nulls , need to remove thm before GWAS


R 
library( "R.utils")
fs <- c( "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR2/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.231005.cpaid.gz", "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR2/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.231005.cpaid.gz", "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR1/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.adjBMI.GenR.EUR.MMG.230731.cpaid.gz","/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR1/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.230912.cpaid.gz") 

## Processing file 
fs <- '/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/GenR1/QCed/ADO/maxLCF/EUR/CLEANED.EGG_cIMT.maxLCF.CLD.MW.noBMI.GenR.EUR.MMG.230912.cpaid.gz'


for (i in 1:length(fs)){

  d <- read.table(gzfile(fs[i]), head=T, string=F)
  file.rename(fs[i], gsub(".cpaid", "_v1.cpaid", fs[i]))
  write.table(d, fs[i], row.names=F, quote=F, col.names=T, sep="\t")
  file.rename(fs[i], gsub(".gz", "", fs[i]))
  gzip(gsub(".gz", "", fs[i]))
}


    
    