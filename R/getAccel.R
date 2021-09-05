
#' Calculate age acceleration.
#'
#' @param c_age
#' A vector includes the chronological age of each sample.
#' @param m_age
#' A vector includes the epigenetic ages with sample order the same with 'c_age'
#' @param method
#' Default: 'Linear', methods to calculate age acceleration, available choices
#' include "None", "Linear", and "Loess".
#' @param do_plot
#' Default: TRUE, whether to visualise the age acceleration results.
#' @param title
#' Default: NA, set the title for the plot.
#' @param point_color
#' Default: NA, define the color for the points in the plot, can be a character or
#' a vector of characters.
#' 
#'
#' @return
#' A vector represents age acceleration in each sample.
#' @export
#' @importFrom
#' stats loess lm setNames cor
#' @examples
#' c_age <- sample(1:100, 20)
#' m_age <- c_age + sample(1:10, 20, replace = TRUE)
#' 
#' ## 1. calculate age acceleration
#' res <- getAccel(c_age, m_age, method='Linear', do_plot=FALSE)
#' print(res)
#' 
#' ## 2. calculate age acceleration and visualise age accel results.
#' pdf('savename.pdf', width=4.3, height=6)
#' res <- getAccel(c_age, m_age, method='Linear', title='TEST')
#' dev.off()
#'
getAccel <- function(c_age, m_age, method='Linear', do_plot=TRUE, title='', point_color=NA){
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
  
  if(do_plot){
    par(mai=c(0.1, 0, 0.05, 0.1))
    plot_legend <- FALSE
    point_pch <- 1
    if(is.na(point_color[1])){
      color <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
      point_pch <- 16
    } else if (length(point_color) > 1){
      color <- as.factor(point_color)
      plot_legend <- TRUE
    } else {
      color <- point_color
    }
    line_col <- 'blue'
    num <- length(c_age)
    if (num < 50){
      cex_size <- 1
    } else if (num < 500){
      cex_size <- 0.8
    } else {
      cex_size <- 0.6
    }
    RMSE <- sqrt(mean((c_age - m_age)^2))
    MAE <- mean(abs(c_age - m_age))
    ## plot m_age vs c_age
    par(fig=c(0, 0.96, 0.45, 0.98), mai=c(0, 0.9, 0.25, 0.1))
    plot(c_age, m_age, xlim=c(min(0, c_age), max(100, c_age)), ylim=c(min(0, m_age), max(100, m_age)), 
         xlab='', ylab='mAge', main=title, cex=cex_size, col=color, pch=point_pch, lwd=1.5)
    
    text(-2, 97, paste0("Pearson's r = ", round(cor(c_age, m_age), 3)),
         cex=0.7, pos=4)
    text(-2, 92,  paste0("RMSE = ", round(RMSE, 2)), cex=0.7, pos=4)
    text(-2, 87,  paste0("MAE = ", round(MAE, 2)), cex=0.7, pos=4)
    if (plot_legend){
      legend(0, 82, levels(color), col=1:length(color), pch=point_pch, box.col='gray', cex=0.6)
    }
    
    ## add fitted line
    if(method == "None"){
      line_col <- 'red'
    }else if(method == "Linear"){
      abline(fit_model, col='red', lty=2, lwd=2)
      text(102, 3, paste0("y = ", round(fit_model$coefficients[2], 2), "x + ",
                         round(fit_model$coefficients[1], 2)), pos=2, col='red', cex=0.7)
    }else if(method == "Loess"){
      j <- order(c_age)
      lines(c_age[j], fit_model$fitted[j], col='red', lty=2, lwd=2)
      text(102, 3, paste0("y = loess_fit(x)"), pos=2, col='red', cex=0.7)
    }
    abline(a=0, b=1, col=line_col, lty=2, lwd=1.5)
    text(100, 8, expression(y==x), pos=2, col=line_col, cex=0.7)
    
    
    ## plot accel vs age
    par(fig=c(0, 0.96, 0, 0.45), mai=c(0.9, 0.9, 0.52, 0.1), new=TRUE)
    plot(c_age, accel, xlim=c(min(0, c_age), max(100, c_age)), 
         ylim=c(min(-30, accel), max(30, accel)), xlab='Chronological Age', 
         ylab='Age Acceleration', cex=cex_size, col=color, pch=point_pch, lwd=1.5)
    abline(a=0, b=0, col='red', lty=2, lwd=1.5)
  }
  
  
  return(accel)
}



