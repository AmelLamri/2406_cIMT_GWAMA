
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/")
  ir_p=paste0(m_p, "/Results/GWAS/CHCP/")
  rw_p=paste0(ir_p, "/Raw/")
  md_p=paste0(ir_p, "/Modified/")
  o_p=paste0(ir_p, "/Cleaned")

# Prepare files 


  for (pheni in "maxRCF"){
    for (sexi in c("MW","WOMEN", "MEN" )){
      for (adji in c("noBMI", "adjBMI") ){
        for(popi in c("EUR")){
          for (agegri in "ADL"){

            filei <- paste0("/CP_cIMT.", pheni, ".", agegri,".",sexi,".",adji,".",popi, ".KL.230922")
             

            if1<- paste0(ir_p,"/Raw/Modified/",pheni, "/",agegri,"/", filei, ".txt")
            if2<- paste0(ir_p,"/Raw/Modified/",pheni, "/",agegri,"/", filei, "_2.txt")
            
            if(!file.exists(paste0(if2, ".gz"))){
              
              if(!file.exists(if1) ){
                unzip( paste0(ir_p,"/Raw/", filei, ".zip"), exdir=paste0(ir_p,"/Raw/Modified/",pheni, "/",agegri,"/") )
              }
              
              cat (" \n reading file \n") 
              chcp <- fread(paste0(if1), head=T)
  
              cat (" \n saving file \n") 
              write.table(chcp,  if2, row.names=F, col.names=T, quote=F, sep="\t" )
              
              cat (" \n zipping file \n") 

              gzip(if2, ext="gz", FUN=gzfile)
              file.remove(if1)

            }

          } 
        }
      }
    }
  }






# Run EasyQC 
  # ADO  
    setwd(paste0(o_p, "/ADO/maxRCF/EUR"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_CHCP_maxRCF_ADO_EUR_QC_quant.ecf"))
     	
  # OAD  EUR
    setwd(paste0(o_p, "/OAD/maxRCF/nonEUR"))
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_OAD_nonEUR_QC_quant.ecf"))

  # OAD  nonEUR
    setwd(paste0(o_p, "/OAD/maxRCF/nonEUR"))
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_OAD_nonEUR_QC_quant.ecf"))

  	
  # ADL  EUR
    setwd(paste0(o_p, "/ADL/EUR/maxRCF/"))
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_ADL_EUR_QC_quant.ecf"))

  # ADL  nonEUR
    setwd(paste0(o_p, "/ADL/maxRCF/nonEUR"))
    EasyQC2(paste0(qc2_p, "//CHCP/EGG_cIMT_CHCP_maxRCF_ADL_nonEUR_QC_quant.ecf"))

