

# m<- "/home/lamria/avr/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results

m_p <- "/home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/"
pheno <- "maxIMT_ADO" 

d <- read.table(paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI.TBL"), head=T, string=F)
nrow(d)
d$chr<-sapply(strsplit(d$MarkerName,":"), `[`, 1)
d$bp<-sapply(strsplit(d$MarkerName,":"), `[`, 2)
d$chr<- as.numeric(d$chr)
d$bp <- as.numeric(d$bp)

jpeg(paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI_manha.jpeg"))
manhattan(d[which(d$chr %in% 1:22),],p="P.value", chr="chr", bp="bp", snp="MarkerName" )
dev.off()


# get SNPnames from TOPMED 

sig <- d[d$P.value< 10^-6, ]


