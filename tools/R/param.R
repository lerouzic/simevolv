
read.param <- function(param.file) {
  # The structure of the parameter file is assumed to be, for each line:
  # PARNAME x y
  # The function returns a list named by the parameter names
  
  lines <- readLines(param.file)
  lines <- lines[!grepl(lines, pattern="^#")]
  splines <- strsplit(lines, split="\\s+")
  nn <- sapply(splines, "[", 1)
  ans <- lapply(splines, function(ll) {
    if (length(ll) <= 1) return(NA)
    ll <- ll[-1]
    numll <- suppressWarnings(as.numeric(ll))
    if (any(is.na(numll))) 
      return(ll)
    else 
      return(numll)
  })
  names(ans) <- nn
  ans
}

write.param <- function(param.file, param) {
  content <- paste(names(param), sapply(param, paste, collapse="\t"), sep="\t")
  writeLines(content, param.file)
}

