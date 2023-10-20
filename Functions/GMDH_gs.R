library(GMDH2)

set.seed(12345)

#########################################################GMDH, model argumanlarında sadece maxlayers ve maxneurons olacak.
GMDH_gs<-function(x, cl, train.percent = 0.75,  maxlayers = 10, maxneurons = 15){
#train.percent : percentage of traning set

x<-t(x)
y <- cl

nobs <- length(y)
indices <- sample(1:nobs)
ntrain <- round(nobs*train.percent,0)
nvalid <- nobs-ntrain

train.indices <- sort(indices[1:ntrain])
valid.indices <- sort(indices[(ntrain+1):nobs])

x.train <- x[train.indices,]
y.train <- y[train.indices]
x.valid <- x[valid.indices,]
y.valid <- y[valid.indices]

# to construct model via GMDH algorithm
model <- GMDH(x.train, y.train, x.valid, y.valid, maxlayers = maxlayers, maxneurons = maxneurons, verbose = FALSE)
return(list(model,model$features))
}
#########################################################modeli ve genleri döndürüyor 

