
#' Calculate age acceleration.
#'
#' @param c_age
#' A vector includes the chronological age of each sample.
#' @param m_age
#' A vector includes the epigenetic ages with sample order the same with 'c_age'
#' @param method
#' Default: 'Linear', methods to calculate age acceleration, available choices
#' include "None", "Linear", and "Loess".
#'
#' @return
#' A vector represents age acceleration in each sample.
#' @export
#' @importFrom
#' stats loess lm setNames
#' @examples
#' c_age <- sample(1:100, 20)
#' m_age <- c_age + sample(1:10, 20, replace = TRUE)
#' res <- getAccel(c_age, m_age, method='Linear')
#' print(res)
#'
getAccel <- function(c_age, m_age, method='Linear'){
  if (length(c_age) < 6){
    method <- "None"
  }
  if (method == "None"){
    message("Age acceleration is calculated as the difference between DNAm age and chronological age.")
    accel <- m_age -c_age
  } else if (method == "Linear") {
    message("Age acceleration is calculated as the residual resulting from a linear regression model which DNAm age is regressed on chronological age.") ## copied
    accel <- lm(m_age ~ c_age)$residuals
  } else if (method == "Loess") {
    message("Age acceleration is calculated as the residual resulting from a nonlinear regression model (loess) which DNAm age is regressed on chronological age. Recommend when sample sie great than 500.")
    accel <- loess(m_age ~ c_age)$residuals
  } else {
    stop(paste0("Method should be one of 'None', 'Linear', and 'Loess'. However got: ", method))
  }
  return(accel)
}
