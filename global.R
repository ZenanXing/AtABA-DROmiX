
display_format <- function(val) {
  if (val != "/" && !is.na(val)) {
    val <- as.numeric(val)
    if (val < 0.01) {
      val <- formatC(val, format = "E", digits = 2)
    } else {
      val <- format(round(val, digits = 2), nsmall = 2)
    }
    return(val)
  } else {
    return("/")
  }
}
