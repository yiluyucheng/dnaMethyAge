######## this script was modified from Daniel Belsky's work (https://github.com/danbelsky/DunedinPACE/blob/main/R/PoAmProjector.R).

preprocessDunedinPACE <- function(betas, ref_means, least_proportion=0.9){
  if (!requireNamespace("preprocessCore", quietly = TRUE)){
    BiocManager::install("preprocessCore")
    require("preprocessCore")
  }else{
    require("preprocessCore")
  }
  common_p <- intersect(rownames(betas), names(ref_means))
  back_p <- length(common_p) / length(ref_means)
  if(back_p > least_proportion){
    betas <- betas[common_p, ]
    betas[,] <- normalize.quantiles.use.target(as.matrix(betas), target=ref_means[common_p])
    missing_p <- setdiff(names(ref_means), common_p)
    #### replace NA with reference mean
    betas <- data.frame(betas)
    betas[missing_p, ] <- NA
    
    ref_means <- ref_means[rownames(betas)]
    for(c in 1:ncol(betas)){
      na_col <- which(is.na(betas[, c]))
      betas[na_col, c] <- ref_means[na_col]
    }
    return(betas)
  }else{
    stop(sprintf("Missing too many probes. Only %.2f%s of the required probes have been found!", back_p * 100, "%"))
  }
}