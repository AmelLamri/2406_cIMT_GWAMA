# DESCIPTION 
#This code complies the differents codes used to run the GWAS meta-analysis of GWAMA 
#path on git /y/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts 
hostname=cat /proc/sys/kernel/hostname
hostname=rigenep1
if[ $(hostname) = rigenep1 ] | [ $(hostname) = rigenep2 ] | [ $(hostname) = ristatp25 ] | [ $(hostname) = ristatp23 ] | [ $(hostname) = ristatp21 ] | [ $(hostname) = ristatp19 ] | [ $(hostname) = ristatp17 ];then 
echo TRUE 
m_p=/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI
metal=/genetics/Programs/Linux-metal/generic-metal/metal 
fi 

o_p=$m_p/Results
s_p=$m_p/Scripts


# Prep files for EasyQC 
  # ABCD

    # in ABCD files, the header line is tab delimited, but the rest of the table is space delimited, which prevents it from being read normally by easyQC
  
    cd /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw
  
    abcdfs=("EGG_cIMT.maxCCA.ADO.MEN.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.maxCCA.ADO.MEN.noBMI.EUR.EM.20230626.txt" "EGG_cIMT.maxCCA.ADO.MW.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.maxCCA.ADO.MW.noBMI.EUR.EM.20230626.txt" "EGG_cIMT.maxCCA.ADO.WOMEN.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.maxCCA.ADO.WOMEN.noBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.MEN.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.MEN.noBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.MW.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.MW.noBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.WOMEN.adjBMI.EUR.EM.20230626.txt" "EGG_cIMT.meanCCA.ADO.WOMEN.noBMI.EUR.EM.20230626.txt")
  
    for filei in ${abcdfs[@]};do 
    
      echo ${abcdfs[@]}
    
      #cp /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/$filei.gz /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei.gz
    
      gunzip /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei.gz
    
      sed -i.bak 's/ /\t/g' /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei
      
      gzip /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei
    
    done
  

  # PANIC
    panicp=/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/PANIC/Raw
    panicpm=$panicp/Modified 
    
    panicf=(PANIC_cIMT.meanLCF.ADO.M.adjBMI.EUR.TS.231213.txt PANIC_cIMT.meanLCF.ADO.M.noBMI.EUR.TS.231213.txt PANIC_cIMT.meanLCF.ADO.MW.adjBMI.EUR.TS.231213.txt PANIC_cIMT.meanLCF.ADO.MW.noBMI.EUR.TS.231218.txt PANIC_cIMT.meanLCF.ADO.W.adjBMI.EUR.TS.231213.txt PANIC_cIMT.meanLCF.ADO.W.noBMI.EUR.TS.231213.txt)
  
    cp $panicp/${panicf[0]}.gz $panicp/tests
    gunzip $panicp/tests/${panicf[0]}.gz
    cat $panicp/tests/${panicf[0]} | head
    #in the study's data, there are quotes, I need to remove them before QC 
    
    for filei in ${panicf[@]};do 
    
      echo ${panicf[@]}
    
      #cp $panicp/$filei.gz $panicpm/$filei.gz
      #cat $panicpm/$filei | head
    
      #gunzip $panicpm/$filei.gz
    
      #sed -i.bak 's/"//g' $panicpm/$filei
      sed -i 's/SNPID/rowN\tSNPID/g' $panicpm/$filei
      #cat $panicpm/$filei | head
    
      gzip $panicpm/$filei
    
    done


  # CHCP 
    chcpf1=(CP_cIMT.maxRCF.ADO.MEN.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.noBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.noBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.noBMI.EUR.KL.230922.zip CP_cIMT.maxRCF.ADO.WOMEN.noBMI.nonEUR.KL.230922)
    
    chcpf2=(CP_cIMT.meanRCF.ADO.MEN.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.noBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.noBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.noBMI.nonEUR.KL.230922)
    
    chcpf12=(CP_cIMT.maxRCF.ADO.MEN.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MEN.noBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.noBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.MW.noBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.adjBMI.EUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.adjBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADO.WOMEN.noBMI.EUR.KL.230922.zip CP_cIMT.maxRCF.ADO.WOMEN.noBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MEN.noBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.MW.noBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.adjBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.adjBMI.nonEUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.noBMI.EUR.KL.230922 CP_cIMT.meanRCF.ADO.WOMEN.noBMI.nonEUR.KL.230922 CP_cIMT.maxRCF.ADL.WOMEN.noBMI.EUR.KL.230922)
    
    chcpp=/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/CHCP/Raw
    chcpm=$chcpp/Modified 
    
    
      #in the study's data, there are quotes, I need to remove them before QC 
      
      for filei in ${chcpf12[@]};do 
      
        echo ${chcpf12[@]}
      
        #cp $chcpp/$filei.zip $chcpm
        #cat $chcpppm/$filei | head
      
         unzip $chcpm/$filei.zip 
        #cat $chcpm/$filei.txt | head
    
        #sed -i.bak 's/"//g' $chcpm/$filei
        #sed -i 's/SNPID/rowN\tSNPID/g' $chcpm/$filei
        #cat $chcpm/$filei | head
      
        gzip $chcpm/$filei.txt
      
      done
  
  # FAMILY 
  




# Run GWAMA for maxIMT 
  cd $o_p/GWAMA 
  $metal $s_p/METAL_MW_maxIMT_NoAdjBMI.metal > $o_p/GWAMA/maxIMT_ADO/METAL_MW_maxIMT_NoAdjBMI.log

# explore GWAMA results 

