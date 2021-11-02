######## this script was modified from Qian Zhang's work(https://github.com/qzhang314/DNAm-based-age-predictor/blob/master/pred.R).


## for each probe, change to missing value to the mean value across all individuals #############
addna <- function(methy){
    mean_p <- mean(methy, na.rm=T)
    if (is.na(mean_p)){
        mean_p <- 0
    }
    methy[is.na(methy)] <- mean_p
    return(methy)
}


preprocessZhangQ2019 <- function(betas){
    sample_id <- rownames(betas)
    if(any(is.na(betas))){
        message("Found NAs, try to replace the NA with mean value for each probe: ")
        betas <- t(apply(betas, 1, addna))
    }
    ### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
    message("Zscore Standardizing:")
    betas <- apply(betas, 2, scale)
    rownames(betas) <- sample_id
    return(betas)
}
