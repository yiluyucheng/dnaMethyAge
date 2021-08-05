######## this script was modified from Qian Zhang's work(https://github.com/qzhang314/DNAm-based-age-predictor/blob/master/pred.R). 



## for each probe, change to missing value to the mean value across all individuals #############
addna <- function(methy){
    methy[is.na(methy)] <- mean(methy,na.rm=T)
    methy[is.na(methy)] <- 0 ##wyc
    return(methy)
}

preprocessZhang2019 <- function(betas){
    message("replace the NA with mean value for each probe ")
    betas <- apply(betas, 2, function(x) addna(x))   

    ### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
    message("Zscore Standardizing:")
    betas <- apply(betas, 1, scale)        
    rownames(betas) <- colnames(betas)
    return(betas)
}




