
require(data.table)
arg <- commandArgs(T)

clock <- arg[1]
beta_file <- arg[2]
info_file <- arg[3]
save_file <- arg[4]
# beta_file <- '/kdmsim1/share/yw19282/projects/EnetTrain/USM_WF_norm_Betas_value.xls'


methyAge <- function(betas, clock='Horvath2013'){
    ## prepare clock coefficients
    usable_clocks <- c('Hannum2013', 'Horvath2013', 'Levine2018', 'Zhang2019')
    if (!(clock %in% usable_clocks)){
        message(paste0("Available clocks are: "))
        message(paste0(usable_clocks, sep=", "))
        stop(paste0("Unavailable for the defined clock: ", clock))
    }
    coefs <- fread(paste0('/home/yw19282/Github/dnaMethAge/coefs/', clock, '.txt'), data.table=FALSE)
    coefs <- setNames(coefs$Coefficient, coefs$Probe)

    ## identify missing probes, set their beta values as zero
    betas <- betas[rownames(betas) %in% names(coefs), ]
    betas[is.na(betas)] <- 0
    missing_probe <- setdiff(names(coefs), c(rownames(betas), "Intercept"))
    if(length(missing_probe) > 0){
        message(paste0("Warnning! Found the below ", length(missing_probe) -1, " probes missing! Will set them to zeros.\n"))
        print(missing_probe)
    }

    ## matrix multiplication
    m_age <- t(betas) %*% matrix(data=coefs[rownames(betas)])
    if("Intercept" %in% names(coefs)){
        m_age <- m_age + coefs["Intercept"]
    }

    if(clock == 'Horvath2013'){
        m_age <- m_age * (20 + 1) + 20
    }

    ## save results
    m_age <- data.frame(Sample=rownames(m_age), mAge=as.numeric(m_age))
    return(m_age)
}


getAccel <- function(c_age, m_age, method='Linear'){
    if (method == "None"){
        message("Age acceleration is computed as the difference between DNAm age and chronological age.")
        accel <- m_age -c_age
    } else if (method == "Linear") {
        message("Age acceleration is computed as the residual resulting from a linear regression model which DNAm age is regressed on chronological age.") ## copied
        accel <- lm(m_age ~ c_age)$residuals
    } else if (method == "Loess") {
        message("Age acceleration is computed as the residual resulting from a nonlinear regression model (loess) which DNAm age is regressed on chronological age. Recommend when sample sie great than 500.") ## copied
        accel <- loess(m_age ~ c_age)$residuals
    } else {
        stop(paste0("Method should be one of 'None', 'Linear', and 'Loess'. However got: ", method))
    }
    return(accel)
}


## read beta file into matrix
betas <- fread(beta_file, data.table=FALSE)
rownames(betas) <- betas$ID_REF
betas <- betas[, -1]

## predict DNAm age based on selected clock
res <- methyAge(betas, clock=clock)

age_accel <- TRUE
# prepare sample information file, should include at least two columns -- Sample, Age.
info <- fread(info_file)
## compute age acceleration
if (age_accel) {
    res <- merge(info, res, by='Sample')
    res$Age_Acceleration <- getAccel(res$Age, res$mAge, method='Linear')
}

fwrite(res, file=save_file, sep='\t')
