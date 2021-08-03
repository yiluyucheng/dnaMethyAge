require(data.table)


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


methyAge <- function(betas, clock='Horvath2013', age_info=FALSE, fit_method='Linear'){
    ## prepare clock coefficients
    usable_clocks <- c('Hannum2013', 'Horvath2013', 'Levine2018')
    if (!(clock %in% usable_clocks)){
        message(paste0("Available clocks are: "))
        message(paste0(usable_clocks, sep=", "))
        stop(paste0("Unavailable for the defined clock: ", clock))
    }
    coefs <- read.table(paste0('../coefs/', clock, '.txt'), header=TRUE)
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
    
    ## post transformation
    if(clock == 'Horvath2013'){
        horvath2013_transform <- function(x){
            if (x > 0){
                x <- x *(20 + 1) + 20
            }else{
                x <- exp(x + log(20 + 1)) - 1
            }        
        }
        m_age[,1] <- sapply(m_age[,1], horvath2013_transform)
    }

    ## save results
    m_age <- data.frame(Sample=rownames(m_age), mAge=as.numeric(m_age))
    
    ## calculate age acceleration
    warning_message <- "\n'age_info' should be a dataframe which contains sample ID and age information, like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge acceleration will not be calculated."
    if (class(age_info) == "logical"){
        if (age_info){
            warning(message(warning_message))
        }
    } else if (class(age_info) == "data.frame") {
            if (all(c('Sample', 'Age') %in% colnames(age_info))){
                m_age <- merge(info, m_age, by='Sample')
                m_age$Age_Acceleration <- getAccel(m_age$Age, m_age$mAge, method=fit_method)
            }else{
                warning(message("\nThe colnames of age_info should include both 'Sample' and 'Age', like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge\nAge acceleration will not be calculated."))
            }
    } else {
        warning(message(warning_message))
    }
    return(m_age)
}


