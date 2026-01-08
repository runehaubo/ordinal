
########################################################
## Change or modify copyright header in R-files:

cp_header <- readLines("~/GitHub/ordinal/ordinal/misc/copyright_header.txt")
cp_src_header <- gsub("#", "/", cp_header)

folder <- "~/GitHub/ordinal/ordinal/R"
filenames <- list.files(folder)

# fn <- filenames[1] # for tests
for(fn in filenames) {
  filepath <- paste(folder, fn, sep="/")
  txt <- readLines(filepath)
  # Get index of copyright header first and last line:
  ind <- grep("^########################################", txt)
  # Check if copyright header exists in file:
  if(grepl("Copyright (c)", txt[ind[1]+1], fixed=TRUE)) {
    txt <- txt[-seq_len(ind[2])] # remove copyright header
    txt <- c(cp_header, txt) # add new copyright header
    writeLines(txt, con=filepath) # write to file
  } else {
    warning(sprintf("No copyright header found in file: %s.", fn))
  }
}
########################################################
## Update header for c-files:

folder <- "~/GitHub/ordinal/ordinal/src"
filenames <- list.files(folder)
# Get *.c and *.h files from /src:
keep <- sapply(strsplit(filenames, ".", fixed = TRUE),
               function(s) s[2] %in% c("c", "h"))
filenames <- filenames[keep]


for(fn in filenames) { # fn <- filenames[1] # for tests
  filepath <- paste(folder, fn, sep="/")
  txt <- readLines(filepath)
  # Get index of copyright header first and last line:
  # ind <- grep("^########################################", txt)
  ind <- grep("^/////////////////////////////////////////", txt)
  # Check if copyright header exists in file:
  if(grepl("Copyright (c)", txt[ind[1]+1], fixed=TRUE)) {
    txt <- txt[-seq_len(ind[2])] # remove copyright header
    txt <- c(cp_src_header, txt) # add new copyright header
    writeLines(txt, con=filepath) # write to file
  } else {
    warning(sprintf("No copyright header found in file: %s.", fn))
  }
}

########################################################
# Write copyright header to new file:
#
# fn <- filenames[1]
for(fn in filenames) {
  filepath <- paste(folder, fn, sep="/")
  txt <- readLines(filepath)
  writeLines(c(cp_header, txt), con=filepath)
  # writeLines(c(cp_src_header, txt), con=filepath)
}

########################################################
