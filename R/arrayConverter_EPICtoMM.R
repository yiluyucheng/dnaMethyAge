#' arrayConverter_EPICtoMM
#' 
#' @description Impute mammalian methylation array data from human Infinium 
#' array data.The coefficients were generated based on a study of n=141 human 
#' blood samples that were profiled using both the mammalian array and the 
#' Illumina 450k array.  each mammalian CpG was fitted by penalized regressions.
#' @param betaEPIC 
#' A beta value dataframe (EPIC or 450k array) with column names as sample ID 
#' and row names as probe name.
#' @param MammalianList 
#' Vector contains probe names from Illumina HorvathMammalianMethylChip40 
#' BeadChip that need to be imputed.
#' @return
#' Imputed beta value matrix.
#' @export
#' @author
#' This is adopted from Lu, Ake T., et al.(PMID: 37563227)
#'
#' @examples
#' data('subGSE174422')
#' n_betas <- arrayConverter_EPICtoMM(betas, c('cg00249943', 'cg00250826', 'cg00292639'))
#' 
arrayConverter_EPICtoMM <- function(betaEPIC, MammalianList){
  message("Converting 450k/EPIC to the mammalian array...")
  data(list='EPICtoMM', envir=environment())
  n_probe <- setdiff(MammalianList, names(fitall))
  if(length(n_probe) > 0){
    stop(paste(c('MammalianList contains', n_probe, 'which are not included in 
                 our training set.'), collapse=" "))
  }
  betaEPIC <- rbind(betaEPIC, Intercept=1)
  ImpMat <- sapply(MammalianList, function(x) .EPICtoMM_one(betaEPIC, x, fitall))
  ImpMat <- t(data.frame(ImpMat, row.names=colnames(betaEPIC)))
  return(ImpMat)
}

.invfn <- function(a){
  2^a/(2^a + 1)
}

.EPICtoMM_one <- function(betaEPIC, probe, fitall){
  coefs <- fitall[[probe]]
  xmat <- betaEPIC[rownames(betaEPIC) %in% coefs$Probe, ,drop=FALSE]
  coefs <- setNames(coefs$Coefficient, coefs$Probe) 
  cg_pred <- t(xmat) %*% matrix(data=coefs[rownames(xmat)])
  return(.invfn(cg_pred))
}


