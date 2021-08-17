### ORIGINAL AUTHOR: Steve Horvath
### Adopted by Yucheng Wang
### modifications: only keep some ensential codes


#source("adjustedBMIQ.R")


imputation <- function(dat1, goldstandard=probeAnnotation21kdatMethUsed$goldstandard2, fastImputation=FALSE){
  #STEP 3: Create the output file called datout
  set.seed(1)
  datMethUsed <- t(dat1[, -1])
  colnames(datMethUsed) <- as.character(dat1[, 1])
  noMissingPerSample <- rowSums(is.na(datMethUsed))
  table(noMissingPerSample)
  max_noMissingPerSample <- max(noMissingPerSample, na.rm=TRUE)
  message("Start data imputation ...")
  #STEP 2: Imputing
  if(! fastImputation & nrow(datMethUsed) > 1 & max_noMissingPerSample < 3000 ){
    # run the following code if there is at least one missing
    if (max_noMissingPerSample > 0){
      if (!requireNamespace("impute", quietly = TRUE)){
        BiocManager::install("impute")
        require("impute")
      }else{
        require("impute")
      }


      dimnames1 <- dimnames(datMethUsed)
      datMethUsed <- data.frame(t(impute.knn(t(datMethUsed))$data))
      dimnames(datMethUsed) <- dimnames1
    } # end of if
  } # end of if (! fastImputation )

  if(max(noMissingPerSample,na.rm=TRUE)>=3000) fastImputation=TRUE

  if(fastImputation | nrow(datMethUsed) == 1){
    if(max_noMissingPerSample >= 3000) {
      normalizeData = FALSE
    }
    # run the following code if there is at least one missing
    if(max_noMissingPerSample >0 & max_noMissingPerSample < 3000){
      dimnames1 <- dimnames(datMethUsed)
      for(i in which(noMissingPerSample>0)){
        selectMissing1 <- is.na(datMethUsed[i,])
        datMethUsed[i, selectMissing1] <- as.numeric(goldstandard[selectMissing1])
      } # end of for loop
      dimnames(datMethUsed) <- dimnames1
    } # end of if
  } # end of if (! fastImputation )
  message("Finished imputation.")
  return(datMethUsed)
}


horvathPreprocess <- function(betas, normalizeData=TRUE){
    #load('horvath_clock.RData')  ## datClock, probeAnnotation21kdatMethUsed
    #probeAnnotation21kdatMethUsed <- read.table('../coefs/27k_reference.txt', header=TRUE)
    data("27k_reference", envir=environment())
    #STEP 2: Restrict the data to 21k probes and ensure they are numeric
    match1 <- match(probeAnnotation21kdatMethUsed$Name , rownames(betas))
    if(sum(is.na(match1)) > 0){
      warning(paste(sum(is.na(match1)), "CpG probes cannot be matched in horvath's ref probes, will set to NA"))
    }
    betas <- merge(probeAnnotation21kdatMethUsed, betas, sort=FALSE, by.x='Name', by.y="row.names", all.x=TRUE)[, -2]
    betas <- imputation(betas, probeAnnotation21kdatMethUsed$goldstandard2)
    # STEP 3: Data normalization (each sample requires about 8 seconds). It would be straightforward to parallelize this operation.
    if(normalizeData){
        message("Normalization by adusted BMIQ ...")
        message(paste0("Estimate running time for normalisation (BMIQ with fixed reference): ", round(5*nrow(betas) / 60, 1), " minutes."))
        betas <- suppressMessages(t(BMIQcalibration(datM=betas, goldstandard.beta=probeAnnotation21kdatMethUsed$goldstandard2, plots=FALSE)))
    }
    return(betas)
}


