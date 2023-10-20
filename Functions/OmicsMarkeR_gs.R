##################### 6) Determan’s Optimal Gene Selection Algorithm

library("OmicsMarkeR")

OmicsMarkeR_gs=function(x, cl, meth=c("rf","svm","glmnet"), f=10,k=3,k.folds=2){ # argumanlar: f, k, k.folds
# possible methods:
#"plsda" (Partial Least Squares Discriminant Analysis), "rf" (Random Forest), "gbm" (Gradient Boosting Machine), "svm" (Support Vector Machines), "glmnet" (Elastic-net Generalized Linear Model), and "pam" (Prediction Analysis of Microarrays)
# plsda, gbm ve pam predictionda problem cikardigi icin almadik

#f : Max.(herhalde) Number of features desired # fs.stability fonksiyonu bu arguman olmadan calismiyor
# k : Number of bootstrapped interations (default=10 uzun surebilir, bu yuzden arguman olarak tanimlayabiliriz)
# k.folds : Number of folds generated during cross-validation (k.folds="LOO" mumkun, leave one out yapiyor)


############### fs.stability: Classification & Feature Selection
fits <- fs.stability(x,
cl,
method = meth,
f = f, # Number of features desired
k = k, # Number of bootstrapped interations
k.folds = k.folds,
verbose = 'none')


result<- new.env()
# frequency (percent identified in all runs).
# we are selecting the ones with frequency = 1:

# Random Forest

ft<-feature.table(fits, meth[1])
result$Features.Selected.by.rf <- as.character(ft$features[ft$freq==1])

# Support Vector Machines

ft2=feature.table(fits, meth[2])
result$Features.Selected.by.svm <- as.character(ft2$features[ft2$freq==1])

# Elastic-net Generalized Linear Model

ft3=feature.table(fits, meth[3])
result$Features.Selected.by.glmnet <- as.character(ft3$features[ft3$freq==1])

res.feature=list(result$Features.Selected.by.rf,result$Features.Selected.by.svm,result$Features.Selected.by.glmnet)

names(res.feature)=c("rf","svm","glmnet")

return(list(fits,res.feature))

}