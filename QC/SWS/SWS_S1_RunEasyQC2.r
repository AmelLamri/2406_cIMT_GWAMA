
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/SWS")
  ir1_p=paste0(m_p, "/Results/GWAS/SWS/")
  rw_p=paste0(ir1_p, "/Raw/")
  md_p=paste0(rw_p, "/Modified/")
  o1_p=paste0(ir1_p, "/QCed/")
 
  # d<- fread("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/SWS/Raw/EGG_cIMT.maxLCF.CLD.MEN.adjBMI.EUR.PT.240130.txt.gz" )


# Prepare files 


  for (pheni in c("maxLCF")){ #
    for (sexi in c("MW", "MEN", "WOMEN")){ # "MW", 
      for (adji in c("noBMI", "adjBMI") ){
        for(popi in c("EUR")){ # , "nonEUR"
          for (agei in c("CLD")){# , "AAD"

            filei <- paste0("/EGG_cIMT.", pheni, ".", agei,".",sexi,".",adji,".",popi, ".PT.240130")
            if1<- paste0(ir1_p,"/Raw/", filei, ".txt")
            of1<- paste0(md_p,"/",  filei, "_2.txt")

            if(!file.exists(paste0(of1, ".gz"))){
              if(!file.exists(paste0(if1, ".gz"))){
                cat( "input files not exist \n ")
  
  
              }else{
                cat (" \n reading file SWS ", pheni, sexi, adji, popi, agei, "\n") 
                sws <- as.data.frame(fread(paste0(if1, ".gz")))
                sws[which(sws$INFO_TYPE=="."),"INFO_TYPE"] <- 0
                #summary(sws[is.na(sws$BETA), "INFO"])
                #summary(sws[which(sws$INFO_TYPE=="0"),"INFO"] )
                cat (" \n            saving file \n") 
                #during QC, there was too many missing betas, so remove them when I save file  
                write.table(sws[!is.na(sws$BETA),], of1, row.names=F, col.names=T, quote=F, sep="\t" )
                cat (" \n            zipping file \n") 
                gzip(of1, ext="gz", FUN=gzfile)
              }
            } 
          }
        }
      }
    }
  }   










 # Run EasyQC 
  # ADO maxLCF EUR 
    setwd(paste0(o1_p, "/ADO/meanLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/SWS_ADO_meanLCF_EUR_QC_quant.ecf"))



