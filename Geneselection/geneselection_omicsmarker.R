
library("plyr")
library("dplyr")
library("caret") 
library("DESeq2")
library("genefilter")
library("foreach")
library("descr")

source("https://raw.githubusercontent.com/mervekasikci/MLGeneSelection/main/Functions/OmicsMarkeR_gs.R")

########################################################################################
################################# Uploading the dataset ################################
########################################################################################

data <- read.table("https://raw.githubusercontent.com/mervekasikci/MLGeneSelection/main/Dataset/KICH_data.txt", header = TRUE)
data <- as.data.frame(data)
data$group <- as.factor(data$group)
data$group <- revalue(data$group, c("1"="NT", "2"="TP"))
y <- filter(data, group == "NT")
x <- y[,-1]


########################################################################################
#################################### Pre-processing ####################################
########################################################################################

############################### Near-zero variance filtering ###########################

out=preProcess(data,method="nzv") 
data1 <- predict(out, data)

################################ Median ratio normalization ############################

data_1<-data.matrix(data1[,-1])
storage.mode(data_1) = "integer"
data_2<-as.matrix(data1[,1])
data_3<-t(data_2)
rownames(data_3)<-"group"
dds <- DESeqDataSetFromMatrix(countData = t(data_1), colData = DataFrame(data_2), ~ as.factor(V1))
dds <- estimateSizeFactors(dds)
normalized_data <- counts(dds, normalized=TRUE)

################################ Logarithmic transformation ############################

transform_data<-log2(normalized_data+1)
nt_data <- rbind(data_3, transform_data)


########################################################################################
######################### Univariate analysis - Student's t-test########################
########################################################################################

x <- nt_data[-1,]
x <- as.matrix(x)
y <- nt_data[1,]
data_new <- cbind(y, t(x))
data_new <- as.data.frame(data_new)
data_new$y <- as.factor(data_new$y)
data_new[,-1] <- sapply(data_new[,-1], as.numeric)
store <- colttests(x=as.matrix(data_new[,-1]), fac=as.factor(data_new[,1]))$p.value
names_genes <- names(data_new)[-1]
data_5 <- cbind(names_genes, store)
data_5 <- as.data.frame(data_5)
data_5$store <- as.numeric(data_5$store)
data_6 <- data_5[order(data_5$store, decreasing = FALSE),] 
data_7 <- data_6[seq(1, 200, 1),]
data_8 <- subset(data_new[,-1], select=data_7$names_genes)
filtered_data <- cbind(y, data_8)


########################################################################################
################################ Development of the model ##############################
########################################################################################

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

data_cv <- filtered_data
size <- dim(data_cv)[1]
k <- 5

set.seed(123) 
indices<-sample(1:size, replace = FALSE)
shuffled_data <- data_cv[indices,]
folds <- cut(seq(1,nrow(shuffled_data)),breaks=k,labels=FALSE)
train.store <- test.store <- list()

for (i in 1:k){
  train.store[[i]] <- shuffled_data[-which(folds==i),]
  test.store[[i]] <- shuffled_data[which(folds==i),]
}

################################### Parallel computing #################################

cl <- makeCluster(5)
registerDoParallel(cl)

result_omicsmarker <- foreach(i = 1:5, .combine=comb, .multicombine = TRUE, 
                              .errorhandling = "remove",
                              .init=list(list(), list(), list()),
                              .packages = c("caret", "OmicsMarkeR", "assertthat")) %dopar% {
                                source("https://raw.githubusercontent.com/mervekasikci/MLGeneSelection/main/Functions/OmicsMarkeR_gs.R")
                                
                                x_train <- train.store[[i]][,-1]
                                x_train <- as.matrix(x_train)
                                storage.mode(x_train) <- "numeric"
                                x_train <- round(x_train, 2)
                                y_train <- train.store[[i]][,1]                                 
                                set.seed(123)
                                model_OmicsMarkeR <- OmicsMarkeR_gs(data.matrix(x_train), as.factor(y_train), meth=c("rf","svm","glmnet"),k=3,k.folds=10)   
                                list(model_OmicsMarkeR[[2]]$rf, model_OmicsMarkeR[[2]]$svm, model_OmicsMarkeR[[2]]$glmnet)
                              } 

getDoParWorkers()                  
stopCluster(cl)

gene_list_rf <- unlist(result_omicsmarker[[1]])
gene_list_rf <- table(gene_list_rf)
gene_list_rf <- as.data.frame(gene_list_rf)
gene_list_rf <- subset(gene_list_rf, Freq>=2) 

gene_list_svm <- unlist(result_omicsmarker[[2]])
gene_list_svm <- table(gene_list_svm)
gene_list_svm <- as.data.frame(gene_list_svm)
gene_list_svm <- subset(gene_list_svm, Freq>=2) 

gene_list_glmnet <- unlist(result_omicsmarker[[3]])
gene_list_glmnet <- table(gene_list_glmnet)
gene_list_glmnet <- as.data.frame(gene_list_glmnet)
gene_list_glmnet <- subset(gene_list_glmnet, Freq>=2) 

write.csv(gene_list_rf, "gene_list_omicsmarker_rf.csv")
write.csv(gene_list_svm, "gene_list_omicsmarker_svm.csv")
write.csv(gene_list_glmnet, "gene_list_omicsmarker_glmnet.csv")
