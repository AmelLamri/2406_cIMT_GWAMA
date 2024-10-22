 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf


  library(EasyQC2)
  library(data.table)
  library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  m0_p="/genetics/Nutrigen/FAMILY/2301_FAMILY_cIMT_GWAS_PHRI"
  ir_p <- paste0(m0_p, "/Outputs")
  qc2_p=paste0(m_p, "/Scripts/QC/FAMILY")
  o_p=paste0(m_p, "/Results/GWAS/FAMILY")



# Prepare files for QC 


  for (pheni in c("maxLCF")){ #
    for (sexi in c("MW","M", "W")){ # "MW", 
      for (adji in c("noBMI", "adjBMI") ){
        for(popi in c("EUR")){ # , "nonEUR"
          for (agegri in c("ADO", "ADL", "OAD", "AAD")){# , "AAD"
            if(agegri != "ADO" & sexi!= "W") next()
            filei <- paste0("/FAMILY_cIMT.", pheni, ".", agegri,".",sexi,".",adji,".",popi, ".AL.20240825")
            if1<- paste0(ir_p,"/", agegri, "/", filei, ".txt")
            if2<- paste0(ir_p,"/", agegri, "/", filei,  "_2.txt")

            if(!file.exists(paste0(if2, ".gz"))){
              
              cat (" \n reading file", pheni, sexi, adji, popi, agegri, "\n") 


              fam <- as.data.frame(fread(paste0(if1, ".gz")))
              colnames(fam)[colnames(fam)=="INFO_TYPE"] <- "INFO_Measure"
              colnames(fam)[colnames(fam)=="IMPUTED"] <- "INFO_TYPE"
              summary(fam[which(fam$INFO_TYPE==0),"INFO"])
              fam[which(fam$INFO_TYPE==0),"INFO"] <- 1
              cat (" \n            saving file \n") 
              write.table(fam, if2, row.names=F, col.names=T, quote=F, sep="\t" )
              cat (" \n            zipping file \n") 
              gzip(if2, ext="gz", FUN=gzfile)
              if(file.exists(if1))file.remove(if1)

            }
          } 
        }
      }
    }
  }





   #ADO 
  setwd(paste0(o_p, "/QCed/ADO/maxLCF/EUR"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_ADO_FAMILY_QC_quant.ecf"))

  #ADL 
  setwd(paste0(o_p, "/QCed/ADL/maxLCF/EUR/"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_ADL_FAMILY_QC_quant.ecf"))

  #OAD
  setwd(paste0(o_p, "/QCed/OAD/maxLCF/EUR"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_OAD_FAMILY_QC_quant.ecf"))
  
  #AAD 
  setwd(paste0(o_p, "/QCed/AAD/maxLCF/EUR"))
  EasyQC2(paste0(qc2_p, "/EGG_cIMT_maxLCF_AAD_FAMILY_QC_quant.ecf"))

