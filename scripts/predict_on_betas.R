
source('MethAge.R')
require(data.table)
arg <- commandArgs(T)

clock <- arg[1]     # Horvath2013 
beta_file <- arg[2]
info_file <- arg[3]
save_file <- arg[4]


if (file.exists(info_file)){
    # prepare sample information file, should include at least two columns -- Sample, Age.
    info <- fread(info_file, data.table=FALSE)    
}else{
    info <- FALSE
}

## read beta file into matrix
# beta_file <- '/kdmsim1/share/yw19282/projects/EnetTrain/USM_WF_norm_Betas_value.xls'
# beta_file <- '/home/yw19282/phd_project/Meth_Age/Horvath_age/MethylationDataExample55.csv'
betas <- fread(beta_file, data.table=FALSE)
rownames(betas) <- betas[, 1]
betas <- betas[, -1]

## predict DNAm age based on selected clock
res <- methyAge(betas, clock=clock, age_info=info, fit_method='Linear', fast_mode=FALSE)

## save results
fwrite(res, file=save_file, sep='\t')
