# Modified from

library(umx)
library(utils)
library(munsell)
library(foreach)
library(doMC)
library(data.table)
library(parallel)
library(R.matlab)

# Resgister your parallel backend. Basically, this means the number in the parantheses should be your number of cores.
registerDoMC(20)

#I recommend these default settings for OpenMx and umx. I find they recover the most accurate estimates.
mxOption(NULL, "Standard Errors", "No")
mxOption(NULL, "Calculate Hessian", "No") 
mxOption(NULL, "Default optimizer", "SLSQP")
umx_set_auto_run(TRUE)
umx_set_auto_plot(FALSE)


# load informations here.
HCP_info <- readMat("../results/genetic_analysis_openmx/input_FG1_FG2_T1wT2w/HCP_OpenMx_info.mat")

# load FG1 FG2 and T1w/T2w
HCP_base_info <- fread("../results/genetic_analysis_openmx/input_FG1_FG2_T1wT2w/hcp_base_info.csv", header=T, data.table=FALSE)

UnivariateACE <- function(D, V, ZYG, MZ, DZ) {
  mzData <- D[D[,ZYG] ==MZ, ]
  dzData <- D[D[,ZYG] ==DZ, ]
  
  m1 <- umxACE(selDVs = V,  sep = "", dzData = dzData, mzData = mzData)
  return(m1)
}

GeneticCorr <- function(A11, A21, A22) {
  res <- NULL
  univariate <- sqrt(A11 * A22)
  bivariate <- A21
  if(univariate == 0) {
    res <- 0
  }
  else {
    res <- bivariate / univariate
  }
  return(res)
}

voxel_num <- 31870
FG1_T1wT2w_HCP <- foreach(i=1:voxel_num, .combine='rbind') %dopar% {
  Analysis <- HCP_base_info
  Analysis$FGone1 <- HCP_info$FGone[, 1, i]
  Analysis$FGone2 <- HCP_info$FGone[, 2, i]
  
  Analysis$T1wT2w1 <- HCP_info$T1wT2w[, 1, i]
  Analysis$T1wT2w2 <- HCP_info$T1wT2w[, 2, i]

  Analysis <- umx_scale(Analysis, varsToScale = c("FGone1", "FGone2", "T1wT2w1", "T1wT2w2"))

  res <- NULL
  res$h1 <- 0
  res$h2 <- 0
  res$h12 <- 0
  res$rG <- 0
  res$Pval <- 1

  if( !(is.nan(mean(Analysis$FGone1)) || is.nan(mean(Analysis$FGone1)) || is.nan(mean(Analysis$T1wT2w1)) || is.nan(mean(Analysis$T1wT2w1)))) {
    Ans_i <- UnivariateACE(Analysis, c("FGone", "T1wT2w"), "Zyg", 1, 2)
    A_mat <- mxEval(top.A, Ans_i)

    res$h1 <- A_mat[1,1]
    res$h2 <- A_mat[2,2]
    res$h12 <- A_mat[2,1]
    res$rG <- GeneticCorr(A_mat[1,1], A_mat[2,1], A_mat[2,2])

    constrained_Ans_i <- umxModify(Ans_i, update = "a_r2c1")
    res_p <-umxCompare(Ans_i, constrained_Ans_i)
    res$Pval <- res_p$p[2]
  }

  return(res)
}
write.csv(FG1_T1wT2w_HCP, file = '../results/genetic_analysis_openmx/output/fg1_t1wt2w_hcp.csv')

FG2_T1wT2w_HCP <- foreach(i=1:voxel_num, .combine='rbind') %dopar% {
  Analysis <- HCP_base_info
  Analysis$FGtwo1 <- HCP_info$FGtwo[, 1, i]
  Analysis$FGtwo2 <- HCP_info$FGtwo[, 2, i]
  
  Analysis$T1wT2w1 <- HCP_info$T1wT2w[, 1, i]
  Analysis$T1wT2w2 <- HCP_info$T1wT2w[, 2, i]
  
  Analysis <- umx_scale(Analysis, varsToScale = c("FGtwo1", "FGtwo2", "T1wT2w1", "T1wT2w2"))
  
  res <- NULL
  res$h1 <- 0
  res$h2 <- 0
  res$h12 <- 0
  res$rG <- 0
  res$Pval <- 1
  
  if( !(is.nan(mean(Analysis$FGtwo1)) || is.nan(mean(Analysis$FGtwo1)) || is.nan(mean(Analysis$T1wT2w1)) || is.nan(mean(Analysis$T1wT2w1)))) {
    Ans_i <- UnivariateACE(Analysis, c("FGtwo", "T1wT2w"), "Zyg", 1, 2)
    A_mat <- mxEval(top.A, Ans_i)
    
    res$h1 <- A_mat[1,1]
    res$h2 <- A_mat[2,2]
    res$h12 <- A_mat[2,1]
    res$rG <- GeneticCorr(A_mat[1,1], A_mat[2,1], A_mat[2,2])
    
    constrained_Ans_i <- umxModify(Ans_i, update = "a_r2c1")
    res_p <-umxCompare(Ans_i, constrained_Ans_i)
    res$Pval <- res_p$p[2]
  }
  
  return(res)
}
write.csv(FG2_T1wT2w_HCP, file = '../results/genetic_analysis_openmx/output/fg2_t1wt2w_hcp.csv')

FG1_FG2_HCP <- foreach(i=1:voxel_num, .combine='rbind') %dopar% {
  Analysis <- HCP_base_info
  Analysis$FGone1 <- HCP_info$FGone[, 1, i]
  Analysis$FGone2 <- HCP_info$FGone[, 2, i]

  Analysis$FGtwo1 <- HCP_info$FGtwo[, 1, i]
  Analysis$FGtwo2 <- HCP_info$FGtwo[, 2, i]

  Analysis <- umx_scale(Analysis, varsToScale = c("FGone1", "FGone2", "FGtwo1", "FGtwo2"))

  res <- NULL
  res$h1 <- 0
  res$h2 <- 0
  res$h12 <- 0
  res$rG <- 0
  res$Pval <- 1
  
  if( !(is.nan(mean(Analysis$FGone1)) || is.nan(mean(Analysis$FGone1)) || is.nan(mean(Analysis$FGtwo1)) || is.nan(mean(Analysis$FGtwo1)))) {
    Ans_i <- UnivariateACE(Analysis, c("FGone", "FGtwo"), "Zyg", 1, 2)
    A_mat <- mxEval(top.A, Ans_i)

    res$h1 <- A_mat[1,1]
    res$h2 <- A_mat[2,2]
    res$h12 <- A_mat[2,1]
    res$rG <- GeneticCorr(A_mat[1,1], A_mat[2,1], A_mat[2,2])    
    constrained_Ans_i <- umxModify(Ans_i, update = "a_r2c1")
    res_p <-umxCompare(Ans_i, constrained_Ans_i)
    res$Pval <- res_p$p[2]
  }

  return(res)
}
write.csv(FG1_FG2_HCP, file = '../results/genetic_analysis_openmx/output/fg1_fg2_hcp.csv')
