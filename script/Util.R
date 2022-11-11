
.mkdir <- function(...) {
    dir.create(..., recursive = TRUE, showWarnings = FALSE)
}

`%&%` <- function(a,b) paste0(a,b)
