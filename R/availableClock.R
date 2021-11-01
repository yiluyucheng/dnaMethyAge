
#' List all available clocks
#'
#' @return A vector includes all available clock names
#' @export
#'
#' @examples 
#' availableClock()
#' 
availableClock <- function(){
  ## list all available clocks
  all_clock <- c('HannumG2013', 'HorvathS2013', 'LevineM2018', 'ZhangQ2019',
                 'ShirebyG2020', 'YangZ2016', 'ZhangY2017')
  all_clock
}