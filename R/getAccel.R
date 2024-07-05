
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
#' @param point_shape
#' Default: NA, define the shape for the points in the plot, can be a single 
#' integer or a vector of characters.
#' @param plot_accel
#' Default: TRUE, whether to plot the accel distribution.
#' @param plot_legend
#' Default: TRUE, whether to plot the legend.
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
getAccel <- function(c_age, m_age, method='Linear', do_plot=TRUE, title='', 
                     point_color=NA, point_shape=NA, plot_accel=TRUE, simple=FALSE,
                     plot_legend=TRUE, x_lim=c(0, 100), y_lim=c(0, 100), x_lab='Chronological Age', 
                     y_lab='mAge', y2_lab='Age Acceleration', marker_cex=0.7, axis_cex=0.9,
                     text_cex=0.7, lab_cex=1, mgp_space=c(3, 1, 0), tck_length=-0.05){
  if (length(c_age) < 6){
    method <- "None"
  }
  if (method == "None"){
    message("Age acceleration is calculated as the difference between DNAm age and chronological age.")
    accel <- m_age -c_age
  } else if (method == "Linear") {
    message("Age acceleration is calculated as the residual resulting from a linear regression model which DNAm age is regressed on chronological age.") ## copied
    fit_model <- lm(m_age ~ c_age, na.action = na.exclude)
    accel <- residuals(fit_model)
  } else if (method == "Loess") {
    message("Age acceleration is calculated as the residual resulting from a nonlinear regression model (loess) which DNAm age is regressed on chronological age. Recommend when sample sie great than 500.")
    fit_model <- loess(m_age ~ c_age, na.action = na.exclude)
    accel <- residuals(fit_model)
  } else {
    stop(paste0("Method should be one of 'None', 'Linear', and 'Loess'. However got: ", method))
  }
  
  if(do_plot){
    
    if(is.na(point_shape[1])){
      point_pch <- 1
    } else if(length(point_shape) == 1 & is.numeric(point_shape)){
      point_pch <- point_shape
    } else {
      point_pch <- as.numeric(as.factor(point_shape))
    }
    base_colors <- colors()
    if(is.na(point_color[1])){
      color <- rgb(red = 0, green = 0, blue = 0, alpha = 0.5)
    } else if (all(point_color %in% base_colors)){
      color <- point_color
    } else {
      color <- as.factor(point_color)
    }
    
    line_col <- 'blue'
    num <- length(c_age)
    
    ylims <- c(min(y_lim[1], m_age), max(y_lim[2], m_age))
    yrange <- ylims[2] - ylims[1]
    xlims <- c(min(x_lim[1], c_age), max(x_lim[2], c_age))
    xrange <- xlims[2] - xlims[1]
    p1_xlab <- x_lab
    
    par(mgp=mgp_space, tck=tck_length)
    if (plot_accel){
      p1_xlab <- ''
      par(mai=c(0.1, 0, 0.05, 0.1))
      par(fig=c(0, 0.96, 0.45, 0.98), mai=c(0, 0.9, 0.25, 0.1))
    }
    
    #### plot m_age vs c_age
    plot(c_age, m_age, xlim=xlims, ylim=ylims, xlab=p1_xlab, ylab=y_lab, main=title, 
          col=color, pch=point_pch, lwd=1.5, cex.lab=lab_cex, cex.axis=axis_cex, cex=marker_cex)

    u_pch <- unique(point_shape)
    if (plot_legend){
      if((length(u_pch) > 1) & (length(levels(color)) > 1)){
        legend(xlims[1], ylims[2] - yrange*0.19, c(levels(color), u_pch), 
               col=c(1:length(levels(color)), rep(1, length(u_pch))), 
               pch=c(rep(16, length(levels(color))), unique(point_pch)), box.col='gray', cex=text_cex)
      }else if(length(levels(color)) > 1){
        legend(xlims[1], ylims[2] - yrange*0.19, levels(color), col=1:length(levels(color)), pch=point_pch[1], box.col='gray', cex=text_cex)
      }else if(length(u_pch) > 1){
        legend(xlims[1], ylims[2] - yrange*0.19, u_pch, col=rep(1, length(u_pch)), pch=unique(point_pch), box.col='gray', cex=text_cex)
      }
    }
    
    ## add fitted line
    if(method == "None"){
      line_col <- 'red'
    }else if(method == "Linear"){
      abline(fit_model, col='red', lty=2, lwd=2)
      text(xlims[2], ylims[1] + yrange*0.03, paste0("y = ", signif(fit_model$coefficients[2], 2), " * x + ",
                                                    signif(fit_model$coefficients[1], 2)), pos=2, col='red', cex=text_cex)
    }else if(method == "Loess"){
      j <- order(c_age)
      lines(c_age[j], fit_model$fitted[j], col='red', lty=2, lwd=2)
      text(xlims[2], ylims[1] + yrange*0.03, paste0("y = loess_fit(x)"), pos=2, col='red', cex=text_cex)
    }
    
    if(!simple){
      text(xlims[1], ylims[2] - yrange*0.03, paste0("Pearson's r = ", signif(cor(c_age, m_age), 3)),
           cex=text_cex, pos=4)
      RMSE <- sqrt(mean((c_age - m_age)^2))
      MAE <- mean(abs(c_age - m_age))
      text(xlims[1], ylims[2] - yrange*0.09,  paste0("RMSD = ", signif(RMSE, 3)), cex=text_cex, pos=4)
      text(xlims[1], ylims[2] - yrange*0.15,  paste0("MAD = ", signif(MAE, 3)), cex=text_cex, pos=4)
      
    }
    abline(a=0, b=1, col=line_col, lty=2, lwd=1.5)
    text(xlims[2], ylims[1] + yrange*0.08, expression(y==x), pos=2, col=line_col, cex=text_cex)
    
    if(plot_accel){
      #### plot accel vs age
      par(fig=c(0, 0.96, 0, 0.45), mai=c(0.9, 0.9, 0.52, 0.1), new=TRUE)
      plot(c_age, accel, xlim=c(min(0, c_age), max(100, c_age)), 
           ylim=c(min(-0.3*yrange, accel), max(0.3 * yrange, accel)), xlab=x_lab, 
           ylab=y2_lab, col=color, pch=point_pch, lwd=1.5, cex.lab=lab_cex, 
           cex.axis=axis_cex, cex=marker_cex)
      
      
      abline(a=0, b=0, col='red', lty=2, lwd=1.5)
    }
  }
  
  
  return(accel)
}



