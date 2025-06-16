
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
                 'ShirebyG2020', 'YangZ2016', 'ZhangY2017', 'LuA2019', 
                 'HorvathS2018', 'DunedinPACE', 'McEwenL2019', 'CBL_specific', 
                 'PCGrimAge', 'PCHorvathS2013', 'PCHannumG2013', 'PCHorvathS2018',
                 'PCPhenoAge', 'CBL_common', 'Cortex_common', 'epiTOC2', 
                 'BernabeuE2023c', 'LuA2023p1', 'LuA2023p2', 'LuA2023p3', 'FuentealbaM2025')
  message("To understand what do these clocks represent for, please refer to 'https://github.com/yiluyucheng/dnaMethyAge'.")
  all_clock
}