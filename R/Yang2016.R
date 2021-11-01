## original author: Yang et al. Genome Biology (2016) 17:205, DOI 10.1186/s13059-016-1064-3
#### EstEpiTOC.R

#### DESCRIPTION
#### An R-function to estimate the pcgtAge-score under the EpiTOC model. Only required input argument is a data matrix of BMIQ (or other type-2 probe correction) normalised DNAm data. Output is the pcgtAge-score for each column (sample) in the input data matrix.

#### INPUT:
#### betas: DNAm data beta-valued matrix with rownames labeling CpGs and columns labeling samples for which a pcgtAge-score is desired
#### mode: this is an optional parameter, which can take on 2 values: "raw" or "z". If the former, function computes the pcgtAge-score as the average methylation over the represented epiTOC loci. If this option is set to "z", then the score is computed as an average of z-scores, reflecting a deviation relative to a specified reference set of samples (as specified by ref.idx)
#### ref.idx: this needs to be specified if mode="z". It is an index vector of column position of samples in betas which should be treated as a reference.


#### OUPTUT:
#### score: a vector of pcgtAge-scores, one value for each sample of betas


EstEpiTOC <- function(betas, epiTOCcpgs, mode="raw", ref.idx=NULL){
  # load("epiTOCcpgs.RData");
  common.v <- intersect(rownames(betas), epiTOCcpgs);
  message(paste0("Number of represented epiTOC CpGs = ", length(common.v), "/", length(epiTOCcpgs)));
  map.idx <- match(common.v, rownames(betas));
  if(mode=="raw"){
    score.v <- colMeans(betas[map.idx,]);
  }
  else if (mode=="z"){
    if(is.null(ref.idx)){
      print("You have selected mode=z, so please specify a reference index vector ref.idx!");
      stop;
      break;
    }
    else {
      sd.v <- apply(betas[map.idx,ref.idx], 1, sd);
      z.m <- (betas[map.idx,] - rowMeans(betas[map.idx,ref.idx])) / sd.v;
      score.v <- colMeans(z.m);
    }
    
  }         
  
  return(score.v);
}