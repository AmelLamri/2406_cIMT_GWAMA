
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



library(EasyQC2)
library(data.table)
library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/CHCP")
  ir_p=paste0(m_p, "/Results/GWAS/CHCP/")
  rw_p=paste0(ir_p, "/Raw/")
  md_p=paste0(ir_p, "/Modified/")
  o_p=paste0(ir_p, "/QCed")

# Prepare files 


  for (pheni in c("maxRCF", "meanRCF", "minRCF",  "maxRCN", "meanRCN", "minRCN")){ #
    for (sexi in c("MW","MEN", "WOMEN")){ # "MW", 
      for (adji in c("noBMI", "adjBMI") ){
        for(popi in c("EUR", "nonEUR")){ # , "nonEUR"
          for (agegri in c("ADO", "ADL", "OAD", "AAD")){# , "AAD"

            filei <- paste0("/CP_cIMT.", pheni, ".", agegri,".",sexi,".",adji,".",popi, ".KL.230922")
            if1<- paste0(ir_p,"/Raw/Modified/", filei, ".txt")
            if2<- paste0(ir_p,"/Raw/Modified/", filei, "_2.txt")
            if3<- paste0(ir_p,"/Raw/Modified/", filei, "_3.txt")

            if(!file.exists(paste0(if3, ".gz"))){
              
              cat (" \n reading file", pheni, sexi, adji, popi, agegri, "\n") 
              if(file.exists(paste0(if2, ".gz"))){

               chcp <- as.data.frame(fread(paste0(if2, ".gz")))
              } else {
                if(!file.exists(if1) ){ 
                  unzip( paste0(ir_p,"/Raw/", filei, ".zip"), exdir=paste0(ir_p,"/Raw/Modified/") )
                  chcp <- as.data.frame(fread(paste0(if1, ".gz")))
                }
                if (file.exists(if1)){
                  chcp <- as.data.frame(fread(paste0(if1)))
                }
              }
              
              chcp[which(chcp$INFO_TYPE==0),"INFO"] <- 1
              cat (" \n            saving file \n") 
              write.table(chcp, if3, row.names=F, col.names=T, quote=F, sep="\t" )
              cat (" \n            zipping file \n") 
              gzip(if3, ext="gz", FUN=gzfile)
              if(file.exists(if1))file.remove(if1)
            }
          } 
        }
      }
    }
  }


tmp<- fread(paste0(if3, ".gz"))

# See who has NA for EAF and why 
table(is.na(tmp$EAF))

summary(tmp$EAF)
tmp[is.na(tmp$EAF),c("N", "N0", "N1", "N2")]



# Run EasyQC 
  # maxRCF 
    # EUR
      # ADO 
        setwd(paste0(o_p, "/ADO/maxRCF/EUR"))
        EasyQC2(paste0(qc2_p, "/CHCP_ADO_maxRCF_EUR_QC_quant.ecf"))
    
      # OAD
        setwd(paste0(o_p, "/OAD/maxRCF/EUR"))
        EasyQC2(paste0(qc2_p, "/CHCP_OAD_maxRCF_EUR_QC_quant.ecf"))
      	
      # ADL
        setwd(paste0(o_p, "/ADL/maxRCF/EUR/"))
        EasyQC2(paste0(qc2_p, "/CHCP_ADL_maxRCF_EUR_QC_quant.ecf"))
    

      # AAD
        setwd(paste0(o_p, "/AAD/maxRCF/EUR/"))
        EasyQC2(paste0(qc2_p, "/CHCP_AAD_maxRCF_EUR_QC_quant.ecf"))
    



    # nonEUR

      # OAD 
        setwd(paste0(o_p, "/OAD/maxRCF/nonEUR"))
        EasyQC2(paste0(qc2_p, "/CHCP_OAD_maxRCF_nonEUR_QC_quant.ecf"))

      # ADL  nonEUR
        setwd(paste0(o_p, "/ADL/maxRCF/nonEUR"))
        EasyQC2(paste0(qc2_p, "//CHCP_ADL_maxRCF_nonEUR_QC_quant.ecf"))  

        # ADO
        setwd(paste0(o_p, "/ADO/maxRCF/nonEUR"))
        EasyQC2(paste0(qc2_p, "/CHCP_ADO_maxRCF_nonEUR_QC_quant.ecf"))

        
      # ADL
        setwd(paste0(o_p, "/AAD/maxRCF/nonEUR/"))
        EasyQC2(paste0(qc2_p, "/CHCP_AAD_maxRCF_nonEUR_QC_quant.ecf"))
    
#HERE< some variants have an EAF of NA because N0==0, N1==0, N2==0, these are found in MW and in MEN samples, they will be repoved during EASYQC

