require(data.table)
arg <- commandArgs(T)

# infile <- '/home/yw19282/phd_project/Meth_Age/deepL/norm_scaled/valid_norm_scaled.xls'
method <- 'None'
method <- arg[1]
infile <- arg[2]
savefile <- arg[3]
 
dat <- fread(infile)

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

dat$Age_Acceleration <- getAccel(dat$Age, dat$mAge)

fwrite(dat, file=savefile, sep=',')

