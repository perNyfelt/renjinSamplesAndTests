
writeFortran <- function(codelist, fname) {
   
   .local <- function(s, chunksize=64) {
      code <- ""
      first <- TRUE
      
      for (i in seq(1, nchar(s), by=chunksize)) {
         code <- paste(code, "     ", ifelse(first, " ", "*"), substr(s, i, i+chunksize-1), "\n", sep="")
         first <- FALSE
      }
      return(code)
   }
   
   code <- paste(lapply(codelist, .local), collapse="")
   
   # Open the file for writing
   fid <- file(fname,"w")
   # Write the code to the file
   writeLines(code, con=fid)
   # Clean up
   close(fid); unlink(fid)
}

WriteOut <- function(md) {
   # Get the current directory
   wd <- getwd()
   # Change to the temporary working directory
   setwd(tempdir())
   # Safety valve - return to the original directory when the function ends.
   on.exit(expr=setwd(wd))
   
   states <- md$states
   inputs <- md$inputs
   outputs <- md$outputs
   pars <- md$pars
   N <- length(states)
   M <- length(inputs)
   S <- length(outputs)
   
   diff.processes <- unique(unlist(lapply(md$sys.eqs, function(x) { names(x$parsed) } )))
   diff.processes <- diff.processes[diff.processes != "dt"]
   
   if (length(diff.processes) != N) { stop( "Number of diffusion processes must match number of states." ) }
   
   sigmat <- vector(mode="list", length=N*N)
   
   k <- 0
   for (i in 1:N) {
      for (j in 1:N) {
         if (is.null(tmp <- md$sys.eqs[[states[j]]]$parsed[[diff.processes[i]]]))
            tmp <- 0
         sigmat[[k<-k+1]] <- tmp
      }
   }
   
   smat <- vector(mode="list", length=S*S)
   
   k <- 0
   for (i in 1:S) {
      for (j in 1:S) {
         k<-k+1
         if (i == j)
            smat[[k]] <- 1
      }
   }
   
   if (length(md$variance) == 0) {
      warning("No variance terms added. Using a diagonal.")
   } else {
      
      for (i in 1:length(md$variance)) {
         if (is.null(tmp <- ParseOutputVariance(md$variance[[i]]$input , outputs))) {
            stop("Failed to understand the provided variance term.")
         }
         smat[[(tmp$i-1) * S + tmp$j]] <- tmp$var
         smat[[(tmp$j-1) * S + tmp$i]] <- tmp$var
      }
   }
   md$obs.eqs$smat <- smat
   
   sigmatvars <- unlist(lapply(unlist(sigmat),all.vars))
   smatvars <- unlist(lapply(unlist(smat),all.vars))
   
   if (md$linear) {
      # Print A, B, C and D
      
      # A
      A <- vector(mode="list", length=N*N)
      
      k <- 0
      for (i in 1:N) {
         for (j in 1:N) {
            A[[k<-k+1]] <- md$sys.eqs[[states[j]]]$A[[states[i]]]
         }
      }
      md$sys.eqs$amat <- A
      
      # B
      B <- vector(mode="list", length=N*M)
      
      if (M > 0) {
         k <- 0
         for (i in seq_len(M)) {
            for (j in seq_len(N)) {
               B[[k<-k+1]] <- md$sys.eqs[[states[j]]]$B[[inputs[i]]]
            }
         }
      }
      md$sys.eqs$bmat <- B
      
      # C
      C <- vector(mode="list", length=S*N)
      
      k <- 0
      for (i in 1:N) {
         for (j in 1:S) {
            C[[k<-k+1]] <- md$obs.eqs[[outputs[j]]]$C[[states[i]]]
         }
      }
      
      md$obs.eqs$cmat <- C
      
      # D
      D <- vector(mode="list", length=S*M)
      
      if (M > 0) {
         k <- 0
         for (i in seq_len(M)) {
            for (j in seq_len(S)) {
               D[[k<-k+1]] <- md$obs.eqs[[outputs[j]]]$D[[inputs[i]]]
            }
         }
      }
      
      md$obs.eqs$dmat <- D
      
#       pars <- unique(c(getPars(A),getPars(B),getPars(C),getPars(D),sigmatvars,smatvars))
#       pars <- unique(c(gsub("(.+)","\\10",states), 
#                        sort(setdiff(pars,c(states,inputs)))))
#       
#       md$NPARAM <- length(pars)
      
      # First generate the A-D matrices.
      code <- c(convertToFortran(A,.amat.aux,N,N,-1,pars,states,inputs),
                convertToFortran(B,.bmat.aux,N,M+1,1,pars,states,inputs),
                convertToFortran(C,.cmat.aux,S,N,-1,pars,states,inputs),
                convertToFortran(D,.dmat.aux,S,M+1,1,pars,states,inputs))
   } else {
      
      # Make the fvec and hvec
      fvec <- vector(mode="list", length=N)
      for (i in 1:N) {
         fvec[[i]] <- md$sys.eqs[[states[i]]]$f
      }
      names(fvec) <- states
      
      hvec <- vector(mode="list", length=S)
      for (i in 1:S) {
         hvec[[i]] <- md$obs.eqs[[outputs[i]]]$h
      }
      names(hvec) <- outputs
      
      fvecxj <- diffVectorFun(fvec,states)
      hvecxj <- diffVectorFun(hvec,states)
      
      if (M != 0) {
         fvecuj <- diffVectorFun(fvec,inputs)
         hvecuj <- diffVectorFun(hvec,inputs)
      } else {
         fvecuj <- NULL
         hvecuj <- NULL
      }
      
#       pars <- unique(c(getPars(fvec),getPars(hvec),sigmatvars, smatvars))
#       pars <- unique(c(gsub("(.+)","\\10",states), 
#                        sort(setdiff(pars,c(states,inputs)))))
#       
#       md$NPARAM <- length(pars)
      
      # First generate the A-D matrices.
      code <- c(convertToFortran(c(fvecxj,fvec),.fvecxj.aux,N,N,0,pars,states,inputs),
                convertToFortran(fvecuj,.fvecuj.aux,N,M+1,1,pars,states,inputs),
                convertToFortran(c(hvecxj,hvec),.hvecxj.aux,S,N,0,pars,states,inputs),
                convertToFortran(hvecuj,.hvecuj.aux,S,M+1,1,pars,states,inputs))
   }
   
   code <- c(code,
             convertToFortran(sigmat,.sigmat.aux,N,N,-1,pars,states,inputs),
             convertToFortran(smat,.smat.aux,S,S,-1,pars,states,inputs))
   
   modelfile <- tempfile(pattern="model")
   writeFortran(code, paste0(modelfile, ".f"))
   writeFortran(genGlobal(N,M,S,md$NPARAM),"global.h")
   
   md$pars <- pars
   
   return( modelfile )
}

ToCovariance <- function(n,data,vnames=NULL) {
   
#    if (n == 1) {
#       warning("Internal error - ToCovariance()  was called with n=1")
#       return(data)
#    }
   
   Nobs <- ifelse(n > 1, nrow(data), length(data))
   cov <- array(dim = c(n,n,Nobs), dimnames=list(vnames, vnames, NULL))
   
   if (n > 1) {
      covtmp <- matrix(nrow=n,ncol=n)
      tril <- row(covtmp) >= col(covtmp)
      seqi <- seq(from=1,to=n*(n+1)/2)
      covtmp[tril] <- seqi
      triu <- match(seqi,t(covtmp))
      
      for (i in 1:Nobs) {
         cov[,,i][tril] <- data[i,]
         cov[,,i][triu] <- data[i,]
      }
   } else {
      cov[1,1,]  <- data
   }
   
   return(cov)
   
}

CheckCorrelationMatrix <- function(corr, vars = rownames(corr)) {
   
   if (!is.matrix(corr))
      stop("Not a matrix.")
   
   # Is corr symmetric?
   if (nrow(corr) != ncol(corr))
      stop("The correlation matrix must be a named square matrix.")
   
   # Named matrix?
   if (is.null(rownames(corr)) || is.null(colnames(corr)))
      stop("The correlation matrix must be a named square matrix.")
   
   if (!(all(vars %in% rownames(corr)) && all(vars %in% colnames(corr))))
      stop("The correlation matrix must contain all MAP variables.")
   
   corr0 <- corr[vars, vars]
   
   # Is the diagonal 1
   if (!all(diag(corr0) == 1))
      stop("The diagonal of the correlation matrix must be 1.")
   
   if (any(abs(corr0) > 1))
      stop("All values in the correlation matrix must between -1 and 1.")
   
   if (!all(corr0 == t(corr0)))
      stop("The correlation matrix is not symmetric.")
   
   tryCatch(chol(corr0), error=function(e) stop("Not a positive semidefinite matrix", call.=FALSE))
   
   return(corr0)
}

IsCompiled <- function(m) {
   if (length(m$libfile) == 0)
      return( FALSE )
   return(file.exists(m$libfile))
}

Compile <- function(m, static = TRUE, verbose = FALSE) {
   # Save the current working directory
   owd = getwd()
   # Revert to the original working directory when exiting
   on.exit(setwd(owd))
   # Change to R's temporary folder
   setwd(wdir <- tempdir())
   
   # Analyse the structure of the model
   m$AnalyseModel()
   # Write out the model Fortran code. The basename (without extension) of the Fortran file is returned: model3242225
   modelfile <- WriteOut(m)
   
   files <- c("covariance.f90", "call_main.c", "mainctsm.f90", "rnorm.c", "simulate.f90", ifelse(m$linear, "dummy_l.f", "dummy_nl.f"))
   
   # Copy the common files
   file.copy(file.path(system.file("common", package = "ctsmr"), files), wdir)
   
   # Delete mainctsm.o if it exists
   if (file.exists("mainctsm.o"))
      unlink("mainctsm.o")
   
   ctsmrdir <- system.file("libctsmr", .Platform$r_arch, package = "ctsmr")
   
   if (static) {
      # Static linking to libRctsm.a
      flag <- normalizePath(file.path(ctsmrdir, "libctsmr.a"))
      if (.Platform$OS.type == "windows") {
         # Convert to the 8+3 format without spaces on Windows. This is a call to a Windows API
         flag <- shortPathName(flag)
         # Place all \\ to /
         flag <- gsub("\\\\", "/", flag)
      }
      flag <- shQuote(flag)
   } else {
      flag <-  sprintf("-L%s -lctsmr -Wl,-rpath,%s", shQuote(ctsmrdir), shQuote(ctsmrdir))
   }
   
   ### Environment variables
   vars <- c("PKG_FFLAGS", "PKG_LIBS", "SED")
   # Save the current value of the environment variables
   envfound <- vars %in% names(Sys.getenv())
   env_saved <- Sys.getenv(vars)
   
   # Set the environment vars
   Sys.setenv(PKG_LIBS = paste0(flag, " $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),
              CYGWIN   = "nodosfilewarning")
   
   libfile <- paste0(modelfile, .Platform$dynlib.ext)
   errorFile <- paste0(basename(libfile), ".err")
   
   # If rtools has been detected, add it to the path only when running R...
   if (!is.null(get_rtools_path())) {
      old <- add_path(get_rtools_path(), 0)
      on.exit(set_path(old), add = TRUE)
   }
   
   cmd <- paste0(shQuote(file.path(R.home(component = "bin"),"/R")),  " CMD SHLIB ",
                 paste(c(files, paste0(basename(modelfile), ".f")), collapse = " "), 
                 " --output=", shQuote(basename(libfile)),
                 " 2> ", shQuote(errorFile))
   
   compiled <- system(cmd, intern = !verbose)
   
   ##### Restore the environment variables
   # Unset those variables that didn't exist.
   Sys.unsetenv(vars[!envfound])
   # Restore the previous value of the already existing variables.
   for (var in vars[envfound]) { eval(paste0("Sys.setenv(", var, "='", env_saved[var], "')")) }
   
   errorMsg <- readLines(errorFile)
   if (verbose)
      writeLines(errorMsg)
   
   if (!file.exists(libfile))
      stop("Compilation error!")
   
   m$libfile <- libfile
   attr(m$libfile, "DLLname") <- basename(modelfile)
}

CallFunction <- function(fun, m, env, verbose = FALSE, unload = FALSE) {
   
   # Check if compiled
   if (!IsCompiled(m)) {
      # Compile the model
      Compile(m, static=TRUE, verbose=verbose)
   }
   
   if (!is.loaded(fun, attr(m$libfile, "DLLname"))) {
      # Load the library
#       lib <- dyn.load(x=m$libfile)
      dyn.load(x=m$libfile)
   }
   
#    fun <- getNativeSymbolInfo(fun, lib)
   ans <- .Call(fun, env, PACKAGE=attr(m$libfile, "DLLname"))
   
   if (unload == TRUE) {
      # Unload again
      dyn.unload(x = m$libfile)
   }
}
