display_format <- function(val) {
  if (is.na(val)) {
    return(NA_character_)
  }
  
  if (val == "/") {
    return("/")
  }
  
  val <- as.numeric(val)
  
  if (val < 0.01) {
    formatC(val, format = "E", digits = 2)
  } else {
    format(round(val, digits = 2), nsmall = 2)
  }
}