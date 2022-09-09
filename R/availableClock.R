
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
                 'ShirebyG2020', 'YangZ2016', 'ZhangY2017', 'Li2018', 'LuA2019', 
                 'HorvathS2018', 'DunedinPACE', 'McEwenL2019', 'CBL_specific', 'CBL_common',
                 'Cortex_common')
  all_clock
}
