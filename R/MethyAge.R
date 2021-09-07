

#' Predict epigenetic age from DNA methylation data
#'
#' @param betas
#' A dataframe with column names as sample ID and row names as probe name,
#' beta = M / (M + U + 100).
#' @param clock
#' Default: 'Horvath2013', define the clock name to use, currently supported
#' clocks are 'Hannum2013', 'Horvath2013', 'Levine2018', and 'Zhang2019'.
#' @param age_info
#' Default: NA, in order to calculate the age acceleration, you need to provide 'age_info'
#' with a dataframe which contains sample ID and age information.
#' @param fit_method
#' Default: 'Linear', select a method to calculate age acceleration, avaliabe
#' choices include "None", "Linear", and "Loess".
#' @param do_plot
#' Default: TRUE, whether to visualise the age acceleration results. Only valid 
#' when age_info is supplied with expected values.
#' @param fast_mode
#' Default: FALSE, whether not to perform data normalisation for the clock of
#' Horvath2013.
#' @param use_cores
#' An integer(e.g. 8) defines the number of cores to use when paralleling preprocess the
#' normalisation step for Horvath2013. Only valid under Linux OS, default value 
#' is the maximum number of cores available.
#'
#' @return
#' A dataframe includes predicated epigenetic ages, age acceleration is only 
#' returned when age_info is supplied with a dataframe contain sample age 
#' information.
#' @export
#' @importFrom
#' utils data
#' @importFrom
#' parallel detectCores mclapply
#'
#'
#' @examples
#' 
#' data('subGSE174422') ## load example dataset, include: betas, info.
#' print(dim(betas))
#' # [1] 485577 8
#' 
#' #### 1. Predict epigenetic age from DNA methylation data
#' horvath_age <- methyAge(betas, clock='Horvath2013', fast_mode=TRUE) ## fast mode
#' 
#' ## Reliable mode, beta values are subjected to fixed-reference based BMIQ 
#' ## normalisation, same as Horvath's paper
#' horvath_age <- methyAge(betas, clock='Horvath2013') 
#' 
#' hannum_age <- methyAge(betas, clock='Hannum2013')
#' 
#' pheno_age <- methyAge(betas, clock='Levine2018')
#' 
#' zhang_age <- methyAge(betas, clock='Zhang2019')
#' 
#' #### 2. Predict epigenetic age and calculate age acceleration
#' hannum_age <- methyAge(betas, clock='Hannum2013', age_info=info)
#' horvath_age <- methyAge(betas, clock='Horvath2013', age_info=info, fit_method='Linear')
#' pheno_age <- methyAge(betas, clock='Levine2018', age_info=info)
#' zhang_age <- methyAge(betas, clock='Zhang2019', age_info=info)
#' 

methyAge <- function(betas, clock='Horvath2013', age_info=NA, fit_method='Linear', 
                     do_plot=TRUE, fast_mode=FALSE, use_cores=detectCores()){
    ## prepare clock coefficients
    usable_clocks <- c('Hannum2013', 'Horvath2013', 'Levine2018', 'Zhang2019')
    if (!(clock %in% usable_clocks)){
        stop(paste(c("Unavailable for the user input clock:", clock,
                     "\n  Please choose one of the available clocks:", usable_clocks), collapse=" "))
    } else if (clock == 'Horvath2013' & !fast_mode){
        #source('preprocessHorvath2013.R')
        betas <- horvathPreprocess(betas, normalizeData=TRUE, use_cores=use_cores)
    } else if (clock == 'Zhang2019'){
        #source('preprocessZhang2019.R')
        betas <- preprocessZhang2019(betas)
    }
    # coefs <- read.table(paste0('../data/', clock, '.txt'), header=TRUE)
    data(list=clock, envir=environment())
    #data('Horvath2013')
    coefs <- setNames(coefs$Coefficient, coefs$Probe)

    ## identify missing probes, set their beta values as zero
    betas <- betas[rownames(betas) %in% names(coefs), ]
    betas[is.na(betas)] <- 0
    missing_probe <- setdiff(names(coefs), c(rownames(betas), "Intercept"))
    if(length(missing_probe) > 0){
        warning(paste(c("Found the below", length(missing_probe),
                        "probes missing! Will set them to zeros.\n ", missing_probe), collapse=" "))
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
    if (class(age_info) == "data.frame") {
            if (all(c('Sample', 'Age') %in% colnames(age_info))){
                m_age <- merge(age_info, m_age, by='Sample', sort=FALSE)
                if (nrow(m_age) < 1){
                    stop(message("Colnames of the input beta dataframe do not match any of the values of the 'Sample' column in age_info!"))
                }
                if("Color" %in% colnames(age_info)){
                    point_color <- m_age$Color
                }else{
                    point_color <- NA
                }
                m_age$Age_Acceleration <- getAccel(m_age$Age, m_age$mAge, method=fit_method, title=clock, do_plot=do_plot, point_color=point_color)
            }else{
                warning(message("\nThe colnames of age_info should include both 'Sample' and 'Age', like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge\nAge acceleration will not be calculated."))
            }
    } else if (is.na(age_info[1])){
        ## age_info is NA, age acceleration will not be calculated.
    } else {
        warning(message(warning_message))
    }
    return(m_age)
}


