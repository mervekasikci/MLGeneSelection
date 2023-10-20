
##################### 4) biosigner
#if (!requireNamespace("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")

#BiocManager::install("biosigner")
library(biosigner)


## signature selection for all 3 classifiers (PLS-DA, RF, SVM)
## we recommend to keep the default bootI = 50 value for your analyzes

biosign_gs<-function(x, cl, bootI=10){ ##arguman olarak sadece bootstrap sayısı alınacak
# bootI:  Number of bootstaps for resampling 
tx<-t(x)
diaSign = biosign(tx, cl, bootI = bootI,pvalN=0.05, fig.pdfC =F,info.txtC="none")
invisible(capture.output(diaSign <- biosign(tx, cl, bootI = bootI,pvalN=0.05, fig.pdfC =F)))
return(list(diaSign,diaSign@"AS"$signatureLs$complete))

}



