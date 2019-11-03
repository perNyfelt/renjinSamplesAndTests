Sys.sleep(5)

## to copy libstan.a 
STATICLIB_EXT <- 'a' 
files <- Sys.glob(paste("*", STATICLIB_EXT, sep = ''))
if (length(files)) { 
   libstanarch <- if (nzchar(R_ARCH)) paste('libctsmr', R_ARCH, sep = '') else 'libctsmr'
   dest <- file.path(R_PACKAGE_DIR, libstanarch)
   message('installing libctsmr.a to ', dest)
   dir.create(dest, recursive = TRUE, showWarnings = FALSE)
   file.copy(files, dest, overwrite = TRUE)
} 