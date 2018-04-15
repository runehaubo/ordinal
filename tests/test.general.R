
txt <- citation("ordinal")
stopifnot(as.logical(grep("year", txt)))
