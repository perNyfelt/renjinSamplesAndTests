#' Validates user given data for CTSM
#' @param data data.frame or list of data.frames.
#' @param inputs vector of names of the input variables.
#' @param outputs vector of names of the output variables.
#' @param firstorderinputinterpolation logical. If TRUE, a first order interpolation between input samples will be used. Otherwise constant hold.
ValidateData <- function(data, inputs, outputs, firstorderinputinterpolation) {
   
   ok <- FALSE
   if (is.data.frame(data)) {
      # Is data a data.frame?
      data <- list(data)
      ok <- TRUE
   } else if (is.list(data)) {
      # Is data a data.frame?
      ok <- all(sapply(data, is.data.frame))
   }
   
   if (!ok)   
      stop("data must be a data.frame or list of data.frames.")
   
   # TODO: This tolerance will be moved to a global package option.
   # Name might be changed too.
   ctsmr.tolerance.constantsampling.detection <- sqrt(.Machine$double.eps)
   
   # Overall sampling time. If all the data.frames are constantly sampled then TRUE.
   # Otherwise FALSE.
   overall.constantsamplingtime <- TRUE
   
   # Loop over all data.frames in the list
   for (i in seq_along(data)) {
      d <- data[[i]]
      # Extract the column names
      cnames <- colnames(d)
      
      # Duplicated variable names?
      if (anyDuplicated(cnames) != 0)
         stop("Duplicated names in data.")
      
      # Is time there? 't' or 'time'
      tname <- cnames[cnames %in% c("t","time")]
      
      # Ambiguous time definition?
      if (length(tname) > 1)
         stop("Both t and time in data. There must only be one time vector.")
      
      t <- d[[tname]]
      
      # Is time increasing?
      if (!identical(t, sort(t)))
         stop(tname, " must be increasing.", domain=NA, call.=FALSE)
      
      # Determine the if the sampling time is constant over the data set
      dt.diff <- diff(t, differences = 2)
      
      constantsamplingtime <- all(abs(dt.diff) <= ctsmr.tolerance.constantsampling.detection)
      
      overall.constantsamplingtime <- overall.constantsamplingtime && constantsamplingtime
      
      ########### INPUT AND OUTPUT CHECKS
      
      AllThere <- function(names, type) {
         idx <- names %in% cnames
         if (any(!idx)) {
            if (sum(!idx) > 1) {
               c1 <- "s "
               c2 <- " are "
            } else {
               c1 <- " "
               c2 <- " is "
            }
            stop(type, c1, paste(names[!idx], collapse=", "), c2, "not present in data.", domain=NA, call.=FALSE)
         }
      }
      
      if (length(inputs) > 0L) {
         ### Check the inputs
         # Are all inputs available?
         AllThere(inputs, "Input")
         
         # NA or 1e300 are not allowed here
         if (any(is.na(d[inputs])) || any(d[inputs] >= 1e300))
            stop("NA are not allowed in the inputs.", call.=FALSE)
      }
      
      ### Check the outputs
      # Are all outputs available?
      AllThere(outputs, "Output")
      
      # Throw a warning if 1e300 if used instead of NA
      if (any(idx <- (d[outputs] >= 1e300), na.rm=TRUE)) {
         warning("Use NA for missing observations rather than a value >= 1e300", domain=NA, call.=FALSE)
         
         # Replace with NA
         d[outputs][idx] <- NA
      }
      
      ########### CHECKS COMPLETE
      
      # Extract the variables we need.
      data[[i]] <- d[c(tname, inputs, outputs)]
      attr(data[[i]], "ConstantSamplingTime") <- constantsamplingtime
      attr(data[[i]], "TimeVar") <- tname
   }
   
   attr(data, "FirstOrderInputInterpolation") <- firstorderinputinterpolation
   attr(data, "ConstantSamplingTime") <- overall.constantsamplingtime
   
   return(data)
}

#' Packs (validated) user given data for CTSM
#' @param data list of data.frames.
#' @param inputs vector of names of the input variables.
#' @param outputs vector of names of the output variables.
PackData <- function(data, inputs, outputs) {
   
   if (!(is.list(data) && all(sapply(data, is.data.frame))))
      stop("data must be a data.frame or list of data.frames.", call.=FALSE, domain=NA)
   
   N <- sapply(data, nrow)
   Nvar <- length(inputs) + length(outputs)
   
   # Allocate the input vector
   input_pack <- numeric(sum(N) * Nvar)
   time_pack <- numeric(sum(N))
   
   time_pack <- unlist(lapply(data, function(d) { d[[1]] }))
   
   input_pack <- unlist(lapply(data, function(d) { c(t(d[-1])) }))
   
   # Replace all NA with 1e300
   input_pack[is.na(input_pack)] <- CTSM_NA_VALUE
   
   attr(input_pack, "time") <- time_pack
   attr(input_pack, "Observations") <- N
   
   return( input_pack )
}
