
#' Imputate array samples, replace NA with sample mean or golden reference
#'
#' @param mt
#' A numeric matrix with column names as sample ID and row names as probe name. 
#' @param ref 
#' Golden reference, a vector consists of numeric values, in which every item 
#' should has a name, e.g. setNames(c(1,2,3), c('a', 'b', 'c'))
#' @param cut_off 
#' Default 0.9, if the ratio of samples with NA in a row greater than this cutoff,
#' then replace all samples with an identical reference value at that row.
#' @param only_ref_rows
#' Default TRUE. TRUE: the final dataframe will only keep rows with names exist 
#' the in reference vector, FALSE: all rows with names presented in original 
#' dataframe and reference vector will be saved.
#' @return
#' A dataframe with all NA replaced.
#' @export
#'
#' @examples
#' 
#' mt <- data.frame(x1=c(1, NA, 3), x2=c(NA, 10, 15), x3=c(3, 5, NA), 
#'                   row.names=c('a', 'b', 'c'))
#' 
#' ref <- ref <- setNames(c(1, 2, 3, 4, 5), c('a', 'b', 'c', 'd', 'e'))
#' 
#' n_mt <- meanImputation(mt, ref)
#' 
#' 
meanImputation <- function(mt, ref, cut_off=0.9, only_ref_rows=TRUE){
  mt <- as.matrix(mt)
  if(only_ref_rows){
    mt <- mt[rownames(mt) %in% names(ref), ]
  }
  t_mt <- rowSums(mt)
  row_mean <- rowMeans(mt[is.na(t_mt),], na.rm=TRUE)
  na_row <- intersect(rownames(mt)[is.na(t_mt)], names(ref))
  for (p in na_row){
    n_miss <- sum(is.na(mt[p,])) 
    if(n_miss / nrow(mt) > cut_off) {
      #mt[p,][is.na(mt[p,])] <- ref[p]
      mt[p,] <- ref[p]
    }else{
      mt[p,][is.na(mt[p,])] <- row_mean[p]
    }
    
  }
  ## replace missing rows with reference values
  miss_row <- setdiff(names(ref), rownames(mt))
  mt <- rbind(mt, matrix(rep(ref[miss_row], ncol(mt)), ncol=ncol(mt), 
                         dimnames = list(miss_row, colnames(mt))))
  
  return(mt)
}