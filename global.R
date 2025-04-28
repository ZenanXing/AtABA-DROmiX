display_format <- function(val) {
  if (val != "/" && !is.na(val)) {
    val <- as.numeric(val)
    return(ifelse(val < 0.01, formatC(val, format = "E", digits = 2), format(round(val, digits = 2), nsmall = 2)))
  } else {
    return(val)
  }
}
