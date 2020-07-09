






parse_functions_source <- function(file_path) {
  
  # Grab list of all functions for export with parallel processing
  is_assign <- function (expr) 
    is.call(expr) && as.character(expr[[1]]) %in% c('=', '<-', 'assign')
  
  is_function <- function (expr) {
    if (! is_assign(expr))
      return(FALSE)
    value = expr[[3]]
    is.call(value) && as.character(value[[1]]) == 'function'
  }
  function_name = function (expr) as.character(expr[[2]])
  function_list <-
    unlist(Map(function_name,   Filter(is_function, parse(file_path))))
  
  return(function_list)
}
# # functions = Filter(is_function, parse("R/row_effect_sizes.R"))
# row_effect_sizes_fun <-
#   unlist(Map(function_name,   Filter(is_function, parse("R/row_effect_sizes.R"))))
# 
# mmd_functions <-
#   unlist(Map(function_name,   Filter(is_function, parse("R/row_effect_sizes.R"))))