#### MASHR STARTING FROM META ANALYSIS ####


library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(zellkonverter)
library(dplyr)
library(edgeR)
library(mashr)


setwd("/sc/arion/projects/CommonMind/mpjanic/PEC/dreamlet_1.3.3/mashr")

library(Matrix)


tabToMatrix <- function(tab, col, rn = "ID", cn = "assay") {
  # check column names
  if (!rn %in% colnames(tab)) stop("Column not found: ", rn)
  if (!cn %in% colnames(tab)) stop("Column not found: ", cn)

  # convert query column names to factor
  if (!is.factor(tab[[rn]])) tab[[rn]] <- factor(tab[[rn]])
  if (!is.factor(tab[[cn]])) tab[[cn]] <- factor(tab[[cn]])

  # extract row and column names for resulting matrix
  rnlvl <- levels(tab[[rn]])
  cnlvl <- levels(tab[[cn]])

  # get indeces
  i <- match(tab[[rn]], rnlvl)
  j <- match(tab[[cn]], cnlvl)

  # convert row,col,value to sparse matrix
  # empty entries are set to 0
  M <- sparseMatrix(i, j,
    x = tab[[col]],
    dims = c(length(rnlvl), length(cnlvl)),
    dimnames = list(rnlvl, cnlvl)
  )

  # convert to real matrix
  data <- as.matrix(M)

  # replace 0's with NA
  # since
  data[data == 0] <- NA

  data
}



run_mash <- function(meta_output) {

tabulator<-read.csv(meta_output)
colnames(tabulator) <- c("ID", "assay", "logFC", "se", "t", "P.Value", "adj.P.Val")

  # convert to matricies
  B <- tabToMatrix(tabulator, "logFC")
  S <- tabToMatrix(tabulator, "se")



  # only keep columns with variance in logFC
  cv <- colVars(B, na.rm = TRUE, useNames = TRUE)
  keep <- (cv > 0) & !is.na(cv)
  B <- B[, keep, drop = FALSE]
  S <- S[, keep, drop = FALSE]

  # run mashr on these matricies
  #-----------------------------

  # set up
  # NA's are replaced with beta = 0 with se = 1e6

B[is.na(B)] <- 0
S[is.na(S)] <- 1e6
  data <- mash_set_data(B, S)

  # estimate model parameters
  U.c <- cov_canonical(data)

  # Estimate correlation structure
  V.em <- mash_estimate_corr_em(data, U.c, details = TRUE)



  # copy model to drop NA terms
  model <- V.em$mash.model

  # B has the same ordering as these, so replace corresponding elements with NA
  # this revents non-sensiccal results for coefficients that were originally NA
  idx <- which(is.na(B))
  model$result$PosteriorMean[idx] <- NA
  model$result$PosteriorSD[idx] <- NA
  model$result$NegativeProb[idx] <- NA
  model$result$lfsr[idx] <- NA

  # format results as new object
  new("dreamlet_mash_result", list(model = model, logFC.original = B))

}


res_mash = run_mash(meta_output='/sc/arion/projects/CommonMind/mpjanic/PEC/dreamlet_1.3.3/final/scz.csv')

# how many genes are significant in at least one cell type
table(get_lfsr(res_mash$model) < 0.05, useNA="ifany")
# how many genes are significant in at least one cell type
table( apply(get_lfsr(res_mash$model), 1, min, na.rm=TRUE) < 0.05)
# how many genes are significant in each cell type
apply(get_lfsr(res_mash$model), 2, function(x) sum(x < 0.05, na.rm=TRUE))


###library(mashr)


#go.gs = get_GeneOntology("CC", to="SYMBOL")
#df_gs = zenith_gsa(res_mash, go.gs)

#write.csv(df_gs, "cerad_mashr_zenith.csv", row.names = FALSE, quote=FALSE)

mashr_output<-get_pm(res_mash$model)

#png(file="cerad_mashr_volcano.png", width=1600, height=1600)
#plotVolcano(res_mash)
#par(mar=c(0,0,0,0))
#dev.off()

#png(file="cerad_mashr_forest_TTR.png", width=600, height=1600)
#plotForest(res_mash, "TTR")
#par(mar=c(0,0,0,0))
#dev.off()

#png(file="cerad_mashr_zenith.png", width=1600, height=1600)
#plotZenithResults(df_gs, 5, 1)
#par(mar=c(0,0,0,0))
#dev.off()

mashr_pval<-get_lfsr(res_mash$model)

write.csv(mashr_output, "scz_mashr_post.csv", quote=FALSE)
write.csv(mashr_pval, "scz_mashr_pval.csv", quote=FALSE)

prob <- compositePosteriorTest(res_mash, colnames(mashr_output))
write.csv(prob, "scz_mashr_comppost.csv", quote=FALSE)
