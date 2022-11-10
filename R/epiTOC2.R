
#' epiTOC2
#' @description Estimate the cumulative number of stem-cell divisions in a sample using the 
#' epiTOC2 model.
#'
#' @param betas 
#' A dataframe with column names as sample ID and row names as probe name,
#' betas = M / (M + U + 100).
#' @param coefs 
#' A dataframe contains coefficents for each probes, colnames as: delta, beta0.
#' @param full_model
#' Default: TRUE, whether to use the full model; In the full model, we use the 
#' estimated methylation ground state values, whereas in the simplified model 
#' we assume that these are all zero. Reason for using the simplified model is 
#' because it could happen that some methylation beta-values in the data matrix 
#' are lower than the ground-state values, which in principle is not allowed. 
#' If this is the case, then more reliable estimates are provided by using the 
#' simplified model.
#' 
#'
#' @return
#' @export
#' @author
#' This function is adopted from Andrew E Teschendorff, 2020 
#' (https://doi.org/10.1186/s13073-020-00752-3) by Yucheng Wang 
#' (wangyucheng511[at]gmail.com)  with modifications.
#'
#' @examples
#' data(list='epiTOC2', package='dnaMethyAge')
#' data('subGSE174422') ## load example dataset, include: betas, info.
#' print(dim(betas))
#' # [1] 485577 8
#' 
#' res <- epiTOC2(betas, coefs, full_model=TRUE)
#' res2 <- epiTOC2(betas, coefs, full_model=FALSE)
#' 
#' 
epiTOC2 <- function(betas, coefs, full_model=TRUE){
  estETOC2.m <- coefs
  data.m <- betas
  ### do epiTOC2
  map.idx <- match(rownames(estETOC2.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""));
  tmp.m <- data.m[map.idx[rep.idx],];
  if(full_model){
    #TNSC.v <- 2 * colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
    TNSC.v <- 2 * colMeans((tmp.m - estETOC2.m[rep.idx,2]) * (1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))),na.rm=TRUE);
  }else{  
    #TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
    TNSC.v <- 2*colMeans(tmp.m * (1/estETOC2.m[rep.idx,1]), na.rm=TRUE);
  }
  return(TNSC.v)
}