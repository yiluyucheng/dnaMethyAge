
require(data.table)
arg <- commandArgs(T)

clock <- arg[1]
beta_file <- arg[2]
save_file <- arg[3]
# beta_file <- '/sharekdmgpu1/share/yw19282/projects/EnetTrain/USM_WF_norm_Betas_value.xls'

## prepare clock coefficients
# clock <- 'Horvath2013'
# clock <- 'Hannum2013'
usable_clocks <- c('Hannum2013', 'Horvath2013', 'Levine2018', 'Zhang2019')
if (!(clock %in% usable_clocks)){
    message(paste0("Available clocks are: "))
    message(paste0(usable_clocks, sep=", "))
    stop(paste0("Unavailable for the defined clock: ", clock))
}
coefs <- fread(paste0('/home/yw19282/Github/dnaMethAge/coefs/', clock, '.txt'), data.table=FALSE)
coefs <- setNames(coefs$Coefficient, coefs$Probe)

## read beta file into matrix
betas <- fread(beta_file, data.table=FALSE)
rownames(betas) <- betas$ID_REF
betas <- betas[, -1]

## identify missing probes, set their beta values as zero
betas <- betas[rownames(betas) %in% names(coefs), ]
betas[is.na(betas)] <- 0
missing_probe <- setdiff(names(coefs), c(rownames(betas), "Intercept"))
if(length(missing_probe) > 0){
    message(paste0("Warnning! Found the below ", length(missing_probe) -1, " probes missing! Will set them to zeros.\n"))
    print(missing_probe)
}

## matrix multiplication
p_age <- t(betas) %*% matrix(data=coefs[rownames(betas)])
if("Intercept" %in% names(coefs)){
    p_age <- p_age + coefs["Intercept"]
}

if(clock == 'Horvath2013'){
    p_age <- p_age * (20 + 1) + 20
}

## save results
p_age <- data.frame(Sample=rownames(p_age), Pred=p_age)
fwrite(p_age, file=save_file, sep=',')

