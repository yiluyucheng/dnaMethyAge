
#' Calculate age acceleration.
#'
#' @param c_age
#' A vector includes the chronological age of each sample.
#' @param m_age
#' A vector includes the epigenetic ages with sample order the same with 'c_age'
#' @param method
#' Default: 'Linear', methods to calculate age acceleration, available choices
#' include "None", "Linear", and "Loess".
#' @param plot_accel
#' Default: TRUE, whether to visualise the age acceleration results.
#' @param title
#' Default: NA, set the title for the plot.
#' 
#'
#' @return
#' A vector represents age acceleration in each sample.
#' @export
#' @importFrom
#' stats loess lm setNames
#' @examples
#' c_age <- sample(1:100, 20)
#' m_age <- c_age + sample(1:10, 20, replace = TRUE)
#' 
#' ## 1. calculate age acceleration
#' res <- getAccel(c_age, m_age, method='Linear', plot_accel=FALSE)
#' print(res)
#' 
#' ## 2. calculate age acceleration and visualise age accel results.
#' pdf('savename.pdf', width=4.3, height=6)
#' res <- getAccel(c_age, m_age, method='Linear', title='TEST')
#' dev.off()
#'
getAccel <- function(c_age, m_age, method='Linear', plot_accel=TRUE, title=''){
  if (length(c_age) < 6){
    method <- "None"
  }
  if (method == "None"){
    message("Age acceleration is calculated as the difference between DNAm age and chronological age.")
    accel <- m_age -c_age
  } else if (method == "Linear") {
    message("Age acceleration is calculated as the residual resulting from a linear regression model which DNAm age is regressed on chronological age.") ## copied
    fit_model <- lm(m_age ~ c_age)
    accel <- fit_model$residuals
  } else if (method == "Loess") {
    message("Age acceleration is calculated as the residual resulting from a nonlinear regression model (loess) which DNAm age is regressed on chronological age. Recommend when sample sie great than 500.")
    fit_model <- loess(m_age ~ c_age)
    accel <- fit_model$residuals
  } else {
    stop(paste0("Method should be one of 'None', 'Linear', and 'Loess'. However got: ", method))
  }
  
  if(plot_accel){
    par(mai=c(0.1, 0, 0.05, 0.1))
    ## plot m_age vs c_age
    par(fig=c(0, 0.96, 0.45, 0.98), mai=c(0, 0.9, 0.25, 0.1))
    color <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
    plot(c_age, m_age, xlim=c(0, 100), ylim=c(0, 100), xlab='', ylab='mAge', 
         main=title, cex=0.6, col=color, pch=16)
    abline(a=0, b=1, col='blue', lty=2, lwd=1.5)
    if(method == "None"){
      abline(a=0, b=1, col='red', lty=2, lwd=2)
    }else if(method == "Linear"){
      abline(fit_model, col='red', lty=2, lwd=2)
    }else if(method == "Loess"){
      j <- order(c_age)
      lines(c_age[j], fit_model$fitted[j], col='red', lty=2, lwd=2)
    }
    
    ## plot accel vs age
    par(fig=c(0, 0.96, 0, 0.45), mai=c(0.9, 0.9, 0.52, 0.1), new=TRUE)
    plot(c_age, accel, xlim=c(0, 100), ylim=c(-30, 30), xlab='Chronological Age', 
         ylab='Age Acceleration', cex=0.6, col=color, pch=16)
    abline(a=0, b=0, col='red', lty=2, lwd=1.5)
  }
  
  
  return(accel)
}



