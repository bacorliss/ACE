



#' Either set working directory to the path of this file, or run this line 
#' within R studio
tryCatch(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))


dir.create(file.path(getwd(), "temp"), showWarnings = FALSE)