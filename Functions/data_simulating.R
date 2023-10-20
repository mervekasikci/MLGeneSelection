library(MASS) #to use rnegbin function. rnegbin generate random outcomes from a Negative Binomial distribution, with mean mu and variance mu + mu^2/theta.
library(edgeR) #to use cpm function. 
data_simulating <- function (number_of_obs = 90, ratio_of_treatment = 0.6, number_of_genes = 5, rate_of_DEG = 0.5, mu=10000, disp=2, fold_change=2, up = 0.5, replace = TRUE) 
{
  m <- number_of_obs/2
  pi0 = 1-rate_of_DEG
  arg <- list(nGenes = number_of_genes, pi0 = pi0, group = c(rep(1, 2* m*ratio_of_treatment), rep(2, 2*m-2*m*ratio_of_treatment)))
  summary(as.factor(arg$group))
  FP <- round(number_of_genes * pi0)
  TP <- number_of_genes - FP
  TP_up <- round(TP * up)
  TP_down <- TP - TP_up
  
  de <- c(rep(1, TP_up), rep(-1, TP_down), rep(0, FP))
  h <- rep(TRUE, number_of_genes)
  counts <- matrix(0, nrow = number_of_genes, ncol = 2 * m)
  delta <- rep(0, number_of_genes)
  if (is.function(fold_change)) {
    lfc <- log(fold_change(TP))
  }
  else {
    lfc <- log(fold_change)
  }
  delta[de != 0] <- lfc * de[de != 0]
  selected_genes <- true_means <- true_disps <- rep(0, number_of_genes)
  left_genes <- 1:length(mu)
  lambda <- phi <- matrix(0, nrow = number_of_genes, ncol = 2 * m)
  while (any(h)) {
    temp <- sample.int(length(left_genes), sum(h), replace)
    temp <- temp[order(temp)]
    selected_genes[h] <- left_genes[temp]
    if (replace == FALSE) {
      left_genes <- left_genes[-temp]
    }
    true_means[h] <- mu[selected_genes[h]]
    true_disps[h] <- disp[selected_genes[h]]
    lambda[h, ] <- matrix(true_means[h], ncol = 1) %*% matrix(rep(1, 
                                                                  2 * m), nrow = 1) * cbind(matrix(rep(1, sum(h) * 
                                                                                                         m), ncol = m), matrix(rep(exp(delta[h]), m), ncol = m))
    phi[h, ] <- matrix(rep(true_disps[h], 2 * m), ncol = 2 * 
                         m)
    counts[h, ] <- rnegbin(sum(h) * 2 * m, lambda[h, ], 1/phi[h, 
    ])
    h <- (rowSums(cpm(counts) > 2) < 3)
  }
  if (any(rowSums(cpm(counts) > 2) < 3)) 
    print("Error: Failed to simulate data: some genes are not expressed.")
  delta <- delta/log(2)
  data <-  t(counts)
  data <- cbind(arg$group, data)
  data <- as.data.frame(data)
  data$V1 <- as.factor(data$V1)
  colnames(data)=c("y",paste0("G",1:(dim(data)[2]-1)))
  return(data = data)
}