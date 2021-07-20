#
# This code is reproduce from Hadley Wickham's devtools package.
# Original: https://raw.github.com/hadley/devtools/master/R/path.r
#
# Used to alter the PATH variable when including Rtools.
#

get_path <- function() {
  strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
}

set_path <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)

  old <- get_path()
  path <- paste(path, collapse = .Platform$path.sep)
  Sys.setenv(PATH = path)
  invisible(old)
}

add_path <- function(path, after = Inf) {
  set_path(append(get_path(), path, after))
}

on_path <- function(...) {
  commands <- c(...)
  stopifnot(is.character(commands))
  unname(Sys.which(commands) != "")
}
