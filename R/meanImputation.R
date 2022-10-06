
#' Imputate array samples, replace NA with sample mean or golden reference
#'
#' @param df
#' A dataframe with column names as sample ID and row names as probe name. 
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
#' df <- data.frame(x1=c(1, NA, 3), x2=c(NA, 10, 15), x3=c(3, 5, NA), 
#'                   row.names=c('a', 'b', 'c'))
#' 
#' ref <- ref <- setNames(c(1, 2, 3, 4, 5), c('a', 'b', 'c', 'd', 'e'))
#' 
#' n_df <- meanImputation(df, ref)
#' 
#' 
meanImputation <- function(df, ref, cut_off=0.9, only_ref_rows=TRUE){
  common_rows <- intersect(rownames(df), names(ref))
  if(only_ref_rows){
    df <- df[rownames(df) %in% names(ref), ]
  }
  for (p in common_rows){
    n_miss <- sum(is.na(df[p,])) 
    if (n_miss > 0){
      if(n_miss / nrow(df) > cut_off) {
        #df[p,][is.na(df[p,])] <- ref[p]
        df[p,] <- ref[p]
      }else{
        df[p,][is.na(df[p,])] <- mean(unlist(df[p,]), na.rm=TRUE)
      }
    }
  }
  ## replace missing rows with reference values
  for (m in setdiff(names(ref), rownames(df))){
    df[m, ] <- ref[m]
  }
  return(df)
}