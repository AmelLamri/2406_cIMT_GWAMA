library(qqman)
library(R.utils)
library(data.table)


# m<- "/home/lamria/avr/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results
#m_p <- "/home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/"#
m_p <- "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI"
rgma <- paste0(m_p,"/Results/GWAMA/")

# Create empty list 
  eumaxp <- list() 
  for (i in 1:4){
    eumaxp[[i]] <- list() 
    for (j in 1:2 ){  
      eumaxp[[i]][[j]] <- list()  
      for (k in 1:3){  
        eumaxp[[i]][[j]][[k]] <- list()  
        }
      }
    }
# Open files tables 

  lft <- read.delim(paste0(m_p, "/Results/FilesTables/241023_GWAS_Files.tbl"), head=T, string=F)
  eumaxpfs <- unique(lft[which(lft$in_ADO_maxLCFp_EUR), c("Study","pop", "pop2", "AgeGroup","sexGroup", "AdjBMI","AdjBMI2","cIMT_measure","QCednoAFO_fn","QCed_p","MetaA_maxp_p"  )])
  nrow(eumaxpfs)
  # lft$AdjBMI2=="noBMI"
  # lft$sexGroup == "MW"
  unique(eumaxpfs[,c("Study","pop2", "AgeGroup","sexGroup", "AdjBMI2")])

fts <- unique( unique(lft[which(lft$in_maxLCFp_EUR), c("pop", "pop2", "AgeGroup","sexGroup", "AdjBMI","AdjBMI2","MetaA_maxp_p"  )]))
fts$ifn <- paste0(fts$MetaA_maxp_p,"/241023_maxLCFp_",fts[,"AgeGroup"],"_",fts[,"pop2"],"_",fts[,"sexGroup"], "_", fts[,"AdjBMI2"] , "_GWAMA_noAFO1.TBL")
fts$ifn2 <- gsub(".TBL","_ChrBp.TBL",fts$ifn)
  agrs <- unique(fts[,"AgeGroup"])
  sgrs <- unique(fts[,"sexGroup"])
  adjs <- unique(fts[,"AdjBMI"])
  pops <- unique(fts[,"pop2"])

 fts[fts$AgeGroup=="ADO","i"]     <- 1
 fts[fts$AgeGroup=="YAD","i"]     <- 2
 fts[fts$AgeGroup=="ADL","i"]     <- 3
 fts[fts$AgeGroup=="OAD","i"]     <- 4
 fts[fts$AgeGroup=="AAD","i"]     <- 5

 fts[fts$AdjBMI2=="noBMI","j"]     <- 1
 fts[fts$AdjBMI2=="adjBMI","j"]     <- 2

 fts[fts$sexGroup=="MW","k"]     <- 1
 fts[fts$sexGroup=="M","k"]     <- 2
 fts[fts$sexGroup=="W","k"]     <- 3

# colnames(fts) <- c(1"AgeGr", 2"SexGr", 3"Adj", 4"abv", 5"i", 6"j",7 "k", 8"ifn", 9"ifn2", 10"N_SNP", 11"N_SNP_nodup")

# Open gwama data and save for locus zoom 
  for (n in 1:nrow(fts)){ #which(fts[,1]=="OAD" 

    if(!file.exists(fts[n,"ifn"])){
      cat ( "file", fts[n,"AgeGroup"], fts[n,"sexGroup"], fts[n,"AdjBMI2"], " not found \n")
      next()
    }else{
      if(class(eumaxp[[fts[n,"i"]]][[fts[n,"j"]]][[fts[n,"k"]]]) == "data.frame") next()
      if(file.exists(paste0(fts[n,"ifn2"], ".gz" ))) {
        cat ( "running 1- ", fts[n,2], fts[n,3], fts[n,4], "maxIMTp \n")
        eumaxp[[fts[n,"i"]]][[fts[n,"j"]]][[fts[n,"k"]]] <- as.data.frame(fread(paste0(fts[n,"ifn2"], ".gz" )))
        tmpord <- eumaxp[[fts[n,"i"]]][[fts[n,"j"]]][[fts[n,"k"]]] 
        names(eumaxp)[fts[n,"i"]] <- fts[n,"AgeGroup"]
        names(eumaxp[[fts[n,"i"]]])[fts[n,"j"]] <-  fts[n,"AdjBMI2"]
        names(eumaxp[[fts[n,"i"]]][[fts[n,"j"]]])[fts[n,"k"]] <- fts[n,"sexGroup"]

 
       } else if(!file.exists(paste0(fts[n,"ifn2"], ".gz" ))){
        cat ( "running 2-", fts[n,"AgeGroup"], fts[n,"sexGroup"], fts[n,"AdjBMI2"], "maxIMTp \n")
        tmpord <- as.data.frame(fread(fts[n,"ifn"]))
        tmpord$Chromosome <- as.numeric(tmpord$Chromosome)
        tmpord$Position<- as.numeric(tmpord$Position)
        tmpord <- tmpord[which((tmpord$Chromosome)%in% 1:22),]
        names(eumaxp)[fts[n,"i"]] <- fts[n,"AgeGroup"]
        names(eumaxp[[fts[n,"i"]]])[fts[n,"j"]] <-  fts[n,"AdjBMI2"]
        names(eumaxp[[fts[n,"i"]]][[fts[n,"j"]]])[fts[n,"k"]] <- fts[n,"sexGroup"]

        
        # re-order chr and position before save for locus zoom
        tmpord$chrpos <- paste0(tmpord$Chromosome, tmpord$Position)
        tmpord <- tmpord[!duplicated(tmpord$chrpos ) & !duplicated(tmpord$chrpos, fromLast=T), ] 
        tmpord <- tmpord[order(tmpord[,"Chromosome"], tmpord[,"Position"]), ] 

        #eumaxp[[fts[n,"i"]]][[fts[n,"j"]]][[fts[n,"k"]]] <- tmpord

        # Save file with chr and position for locus zoom 

        write.table(tmpord, fts[n,"ifn2"] , row.names=F, col.names=T, quote=F, sep="\t")
        #if(!file.exists(paste0(fts[n,"ifn2"], ".gz"))) 
        gzip (fts[n,"ifn2"])
      }
    }
  }

# Plot Manhattan 

      tmpord$"P-value" <- as.numeric(tmpord$"P-value")
      if(!file.exists(gsub("_ChrBp.TBL.gz","_Manh.jepg", fts[n,"ifn2"]))){
        jpeg(gsub("_ChrBp.TBL.gz","_Manh.jepg", fts[n,"ifn2"]))
          manhattan(tmpord,p="P-value", chr="Chromosome", bp="Position", snp="MarkerName" )
        dev.off()
 
      }
`
# Keep only significant loci 
  sig <- eumaxp

  for (i in 1:4){
    for (j in 1 ){  
      for (k in 1:2){  
        cat (names(eumaxp)[i] , names(eumaxp[[i]])[j] , names(eumaxp[[i]][[j]])[k] , "\n" ) 
        colnames(sig[[i]][[j]][[k]]) <- gsub("P.value","P-value", colnames(sig[[i]][[j]][[k]]))
        sig[[i]][[j]][[k]] <- sig[[i]][[j]][[k]][sig[[i]][[j]][[k]][,"P-value",] <= 0.05,]
        sig[[i]][[j]][[k]] <- sig[[i]][[j]][[k]][order(sig[[i]][[j]][[k]][,"P-value"]),]

      }
    }
  }

  sig2 <- sig
   for (i in 1:4){
      for (j in 1 ){  
        for (k in 1:2){  
          colnames(sig2[[i]][[j]][[k]]) [! colnames(sig2[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")] <- paste0(colnames(sig2[[i]][[j]][[k]]) [! colnames(sig2[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")], "_",names(sig2)[i],"_",names(sig2[[i]])[[j]],"_",names (sig2[[i]][[j]])[[k]])
          sig2[[i]][[j]][[k]] <- sig2[[i]][[j]][[k]][order(sig2[[i]][[j]][[k]]$"P-value" ),]
        }
      }
   }
  
  head(sig2[[3]][[1]][[1]])
  head(sig2[[4]][[1]][[1]])
  
  min(sig2[[1]][[1]][[1]]$"P-value" )


summary(sig[[4]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    491     491     864    1295    2114    3469
summary(sig[[3]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    657    1779    2114    2395    2771    4550

summary(sig[[2]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    302    1575    1575    1383    1877    1877
summary(sig[[1]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    241    1650    2752    2729    3875    5284




# EUR MW NoBMI 



# Keep only significant loci 
  eumaxpa <- eumaxp

  for (i in 1:4){
    for (j in 1 ){  
      for (k in 1:2){  
        cat (names(eumaxpa)[i] , names(eumaxpa[[i]])[j] , names(eumaxpa[[i]][[j]])[k] , "\n" ) 
        colnames(eumaxpa[[i]][[j]][[k]]) <- gsub("P.value","P-value", colnames(sig[[i]][[j]][[k]]))
          colnames(eumaxpa[[i]][[j]][[k]]) [! colnames(eumaxpa[[i]][[j]][[k]]) %in% c("MarkerName")] <- paste0(colnames(eumaxpa[[i]][[j]][[k]]) [! colnames(eumaxpa[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")], "_",names(eumaxpa)[i],"_",names(eumaxpa[[i]])[[j]],"_",names (eumaxpa[[i]][[j]])[[k]])


      }
    }
  }


pci <- c("P-value_ADO_MW_noBMI","P-value_YAD_MW_noBMI","P-value_ADL_MW_noBMI","P-value_OAD_MW_noBMI")
eci <- c("Effect_ADO_MW_noBMI","Effect_YAD_MW_noBMI", "Effect_ADL_MW_noBMI","Effect_OAD_MW_noBMI")
epci <-  c("MarkerName","Effect_ADO_MW_noBMI","P-value_ADO_MW_noBMI","Effect_YAD_MW_noBMI","P-value_YAD_MW_noBMI", "Effect_ADL_MW_noBMI","P-value_ADL_MW_noBMI","Effect_OAD_MW_noBMI","P-value_OAD_MW_noBMI")



  MW_noBMI_2 <- merge(eumaxpa[[1]][[1]][[1]], eumaxpa[[2]][[1]][[1]], all=T)
  MW_noBMI_2 <- merge(MW_noBMI_2           , eumaxpa[[3]][[1]][[1]], all=T)
  MW_noBMI_2 <- merge(MW_noBMI_2           , eumaxpa[[4]][[1]][[1]], all=T)

MW_noBMI_2a <- MW_noBMI_2[order(MW_noBMI_2$"P-value_OAD_MW_noBMI"),]
write.csv(MW_noBMI_2a[MW_noBMI_2a$"P-value_OAD_MW_noBMI"<0.0005,epci], paste0(rgma, "/OAD_Sig_loci_across_Ages.csv"), row.names=F, quote=F)


MW_noBMI_2b<- MW_noBMI_2[order(MW_noBMI_2$"P-value_ADL_MW_noBMI"),]
write.csv(MW_noBMI_2b[which(MW_noBMI_2$"P-value_ADL_MW_noBMI"<0.0005),epci], paste0(rgma, "/ADL_Sig_loci_across_Ages.csv"), row.names=F, quote=F)


MW_noBMI_2c<- MW_noBMI_2[order(MW_noBMI_2$"P-value_YAD_MW_noBMI"),]
write.csv(MW_noBMI_2c[which(MW_noBMI_2$"P-value_YAD_MW_noBMI"<0.005),epci], paste0(rgma, "/YAD_Sig_loci_across_Ages.csv"), row.names=F, quote=F)

MW_noBMI_2d<- MW_noBMI_2[order(MW_noBMI_2$"P-value_ADO_MW_noBMI"),]
write.csv(MW_noBMI_2d[which(MW_noBMI_2$"P-value_ADO_MW_noBMI"<0.005),epci], paste0(rgma, "/ADO_Sig_loci_across_Ages.csv"), row.names=F, quote=F)




colnames(MW_noBMI_2) <- gsub("_ADO_MW_noBMI_ADO_MW_noBMI", "_ADO_MW_noBMI", colnames(MW_noBMI_2))
colnames(MW_noBMI_2) <- gsub("_ADL_MW_noBMI_ADL_MW_noBMI", "_ADL_MW_noBMI", colnames(MW_noBMI_2))

colnames(MW_noBMI_2) <- gsub("_OAD_MW_noBMI_OAD_MW_noBMI", "_OAD_MW_noBMI", colnames(MW_noBMI_2))

colnames(MW_noBMI_2) <- gsub("_YAD_MW_noBMI_YAD_MW_noBMI", "_YAD_MW_noBMI", colnames(MW_noBMI_2))

