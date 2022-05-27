

#' Predict epigenetic age from DNA methylation data
#'
#' @param betas
#' A dataframe with column names as sample ID and row names as probe name,
#' beta = M / (M + U + 100).
#' @param clock
#' Default: 'HorvathS2013', define the clock name to use, currently supported
#' clocks are 'HannumG2013', 'HorvathS2013', 'LevineM2018', and 'ZhangQ2019'.
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
#' HorvathS2013.
#' @param use_cores
#' An integer(e.g. 8) defines the number of cores to use when paralleling preprocess the
#' normalisation step for HorvathS2013. Only valid under Linux OS, default value 
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
#' horvath_age <- methyAge(betas, clock='HorvathS2013', fast_mode=TRUE) ## fast mode
#' 
#' ## Reliable mode, beta values are subjected to fixed-reference based BMIQ 
#' ## normalisation, same as Horvath's paper
#' horvath_age <- methyAge(betas, clock='HorvathS2013') 
#' 
#' hannum_age <- methyAge(betas, clock='HannumG2013')
#' 
#' pheno_age <- methyAge(betas, clock='LevineM2018')
#' 
#' zhang_age <- methyAge(betas, clock='ZhangQ2019')
#' 
#' #### 2. Predict epigenetic age and calculate age acceleration
#' hannum_age <- methyAge(betas, clock='HannumG2013', age_info=info)
#' horvath_age <- methyAge(betas, clock='HorvathS2013', age_info=info, fit_method='Linear')
#' pheno_age <- methyAge(betas, clock='LevineM2018', age_info=info)
#' zhang_age <- methyAge(betas, clock='ZhangQ2019', age_info=info)
#' 



methyAge <- function(betas, clock='HorvathS2013', age_info=NA, fit_method='Linear', 
                     do_plot=TRUE, fast_mode=FALSE, use_cores=detectCores()){
    ## prepare clock coefficients
    usable_clocks <- availableClock()
    if (!(clock %in% usable_clocks)){
        stop(paste(c("Unavailable for the user input clock:", clock,
                     "\n  Please choose one of the available clocks:", usable_clocks), collapse=" "))
    } else {
        # coefs <- read.table(paste0('../data/', clock, '.txt'), header=TRUE)
        data(list=clock, envir=environment())
        is_beta <- TRUE
        plot_simple <- FALSE
        x_lim = y_lim = c(0, 100)
        if (clock == 'HorvathS2013' & !fast_mode){
            betas <- horvathPreprocess(betas, normalizeData=TRUE, use_cores=use_cores)
        } else if (clock == 'ZhangQ2019'){
            betas <- preprocessZhangQ2019(betas)
        } else if (clock == 'YangZ2016'){
            mAge <- EstEpiTOC(betas, epiTOCcpgs, mode="raw", ref.idx=NULL)
            m_age <- as.matrix(mAge)
            is_beta <- FALSE
        } else if(clock == 'DunedinPACE'){
            betas <- preprocessDunedinPACE(betas, ref_means=gold_standard_means)
        }
        ## Free the Y limits in plotting
        if(clock %in% c('YangZ2016', 'ZhangY2017', 'LuA2019')){
            y_lim = c(NULL, NULL)
            plot_simple <- TRUE
        }
        
    }
    
    if(is_beta){
        #data('HorvathS2013')
        coefs <- setNames(coefs$Coefficient, coefs$Probe)
        ## add intercept
        betas <- rbind(betas, Intercept=rep(1, ncol(betas)))
        
        ## identify missing probes, set their beta values as zero
        betas <- betas[rownames(betas) %in% names(coefs), ]
        betas[is.na(betas)] <- 0
        missing_probe <- setdiff(names(coefs), rownames(betas))
        if(length(missing_probe) > 0){
            warning(paste(c("Found the below", length(missing_probe),
                            "probes missing! Will set them to zeros.\n ", missing_probe), collapse=" "))
        }
    
        ## matrix multiplication
        m_age <- t(betas) %*% matrix(data=coefs[rownames(betas)])
    
        ## post transformation
        if(clock %in% c('HorvathS2013', 'ShirebyG2020', 'HorvathS2018', 'McEwenL2019')){
            HorvathS2013_transform <- function(x){
                if (x > 0){
                    x <- x *(20 + 1) + 20
                }else{
                    x <- exp(x + log(20 + 1)) - 1
                }
            }
            m_age[,1] <- sapply(m_age[,1], HorvathS2013_transform)
        } else if(clock %in% c('CBL_specific', 'CBL_common', 'Cortex_common')){
            m_age[,1] <- exp(m_age[,1])
        }
        
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
                if("Shape" %in% colnames(age_info)){
                    point_pch <- m_age$Shape
                }else{
                    point_pch <- NA
                }
                m_age$Age_Acceleration <- getAccel(m_age$Age, m_age$mAge, method=fit_method, title=clock, do_plot=do_plot, point_color=point_color, point_shape=point_pch, simple=plot_simple, x_lim=x_lim, y_lim=y_lim)
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




