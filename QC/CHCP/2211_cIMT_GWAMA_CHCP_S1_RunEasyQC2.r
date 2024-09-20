
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



  library(EasyQC2)
  library(data.table)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/")
  ir_p=paste0(m_p, "/Results/GWAS/CHCP/")
  rw_p=paste0(ir_p, "/Raw/")
  md_p=paste0(ir_p, "/Modified/")
  o_p=paste0(ir_p, "/Results/GWAS/CHCP/Cleaned")

# Prepare files 


  files <- c("", "CP_cIMT.maxRCF.OAD.MEN.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.noBMI.EUR.KL.230922")

  for(pheni in "maxRCF"){
  for(agri in "maxRCF"){

  for(i in 1:length(files)){
    cat (i, "\n")
    filei <- paste0("/CP_cIMT.", pheni, ".OAD.MEN.adjBMI.nonEUR.KL.230922")
    unzip( paste0(r_p,"/", filei, ".zip"), exdir=paste0(r_p,"/Modified") )
    chcp <-  fread(paste0(r_p,"/Modified/",pheni, "/", filei, ".txt"), head=T)
    write.table(chcp,paste0(r_p,"/Modified/",pheni, "/",filei, ".txt"), row.names=F, col.names=T, quote=F, sep="\t" )
    gzip(paste0(r_p,"/Modified/", pheni, "/", filei, ".txt"), ext="gz", FUN=gzfile)

    }





# Run EasyQC 
  # ADO  
    setwd(paste0(o_p, "/ADO/maxRCF/EUR"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_CHCP_maxRCF_ADO_EUR_QC_quant.ecf"))
     	
  # OAD  
    setwd(paste0(o_p, "/OAD/maxRCF/EUR"))
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_OAD_EUR_QC_quant.ecf"))

