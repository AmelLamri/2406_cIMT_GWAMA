

r_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw"
#easy QC2 does not recognize .zip files, and cannot open the txt files as they are (maybe because space delimited ? )



#OAD 
for fi in CP_cIMT.maxRCF.OAD.MEN.adjBMI.EUR.KL.230922 /CP_cIMT.maxRCF.OAD.MEN.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.OAD.MW.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.OAD.WOMEN.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.OAD.WOMEN.noBMI.EUR.KL.230922;do
 unzip $r_p/$fi.zip -d $r_p
 done




R 
r_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw"
library(data.table)

# OAD EUR
files <- c("CP_cIMT.maxRCF.OAD.MEN.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MEN.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.noBMI.EUR.KL.230922")
for(i in 1:length(files)){
  chcp <-  fread(paste0(r_p,"/", files[i], ".txt"), head=T)
  write.table(chcp,paste0(r_p,"/", files[i], ".txt"), row.names=F, col.names=T, quote=F, sep="\t" )
  gzip(paste0(r_p,"/", files[i], ".txt"), ext="gz", FUN=gzfile)
}





  # OAD nonEUR
    r_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw"

    files <- c("CP_cIMT.maxRCF.OAD.MEN.adjBMI.nonEUR.KL.230922", "CP_cIMT.maxRCF.OAD.MEN.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.noBMI.EUR.KL.230922")


    for(i in 1:length(files)){
      cat (i, "\n")
      unzip( paste0(r_p,"/", files[i], ".zip"), exdir=paste0(r_p,"/Modified") )
      chcp <-  fread(paste0(r_p,"/Modified/", files[i], ".txt"), head=T)
      write.table(chcp,paste0(r_p,"/Modified/", files[i], ".txt"), row.names=F, col.names=T, quote=F, sep="\t" )
      gzip(paste0(r_p,"/Modified/", files[i], ".txt"), ext="gz", FUN=gzfile)

    }


  # ADO nonEUR
    r_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw"

    files <- c("CP_cIMT.maxRCF.OAD.MEN.adjBMI.nonEUR.KL.230922", "CP_cIMT.maxRCF.OAD.MEN.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.MW.noBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.adjBMI.EUR.KL.230922", "CP_cIMT.maxRCF.OAD.WOMEN.noBMI.EUR.KL.230922")

for (agrgri in c("ADO", "OAD","ADL", "AAD")){
    for(i in 1:length(files)){
      cat (i, "\n")
      unzip( paste0(r_p,"/", files[i], ".zip"), exdir=paste0(r_p,"/Modified") )
      chcp <-  fread(paste0(r_p,"/Modified//maxRCF/", agrgri", files[i], ".txt"), head=T)
      write.table(chcp,paste0(r_p,"/Modified/ADO/", files[i], ".txt"), row.names=F, col.names=T, quote=F, sep="\t" )
      gzip(paste0(r_p,"/Modified/ADO/", files[i], ".txt"), ext="gz", FUN=gzfile)

    }

}
