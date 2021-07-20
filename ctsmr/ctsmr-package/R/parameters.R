
ctsm$methods(
   
   setParameter = function(...) {
      'Use to set the initial value and lower and upper bounds.'
      
      if (nargs() == 0L) { stop("You must give me some names vectors to work with.") }
      
      cc <- list(...)
      
      var.names <- names(cc)
      
      # All elements in cc should be named
      if (!all(nzchar(var.names))) { stop("All arguments must be of the form f=c(initial=).") }
      
      if (!all(table(var.names) == 1L)) { stop("Multiple definitions of the same parameter.") }
      
      if (!all(sapply(cc, is.numeric))) { stop("All arguments must be numeric vectors.") }
      
      for (var in var.names) {
         
         # Need to check what we know about var
         
         gotInitial <- FALSE
         gotLower <- FALSE
         gotUpper <- FALSE
         gotPrior <- FALSE
         
         lower <- upper <- prior <- NA
         
         par <- cc[[var]]
         par.names <- names(par)
         
         par.named <- par[i <- nzchar(par.names)]
         if (is.null(par.names))
            par.unnamed <- par
         else
            par.unnamed <- par[!i]
         
         ok <- TRUE
         
         for (name in names(par.named)) {
            
            ok <- TRUE
            
            if (name %in% c("init","initial")) {
               if (gotInitial) { 
                  warning(sprintf("Multiple instances of the initial value for %s", var))
                  ok <- FALSE
               } else {
                  initial <- par.named[[name]]
                  gotInitial <- TRUE
               }
               
            } else if (name %in% c("lower","lowerbound","lb")) {
               if (gotLower) { 
                  warning(sprintf("Multiple instances of the lowerbound value for %s", var))
                  ok <- FALSE
               } else {
                  lower <- par.named[[name]]
                  gotLower <- TRUE
               }
            } else if (name %in% c("upper","upperbound","ub")) {
               if (gotUpper) { 
                  warning(sprintf("Multiple instances of the upperbound value for %s", var))
                  ok <- FALSE
               } else {
                  upper <- par.named[[name]]
                  gotUpper <- TRUE
               }
            } else if (name %in% c("prior","priorsd","psd")) {
               if (gotPrior) { 
                  warning(sprintf("Multiple instances of the prior value for %s", var))
                  ok <- FALSE
               } else {
                  prior <- par.named[[name]]
                  gotPrior <- TRUE
               }
            } else {
               warning(sprintf("Don't know what to do with %s for %s", name, var))
               ok <- FALSE
            }
            if (!ok) { break }
         }
         if (!ok) { break }
         
         n <- length(par.unnamed)
         i <- 0
         
         if (!gotInitial && i < n) {
            initial <- par.unnamed[i <- i + 1]
            gotInitial <- TRUE
         }
         
         if (!gotLower && i < n) {
            lower <- par.unnamed[i <- i + 1]
            gotLower <- TRUE
         }
         
         if (!gotUpper && i < n) {
            upper <- par.unnamed[i <- i + 1]
            gotUpper <- TRUE
         }
         
         if (!gotPrior && i < n) {
            prior <- par.unnamed[i <- i + 1]
            gotPrior <- TRUE
         }
         
         if (n > i) { warning("The additional value (", paste(par.unnamed[-seq(1,i)],collapse=","), ") given for ", var, " are ignored.") }
         
         if (!gotInitial) {
            warning("No initial value for ", var)
            next
         }
         
         if (gotLower + gotUpper == 1) {
            warning("Cannot find the ", ifelse(gotLower, "upperbound", "lowerbound"), " for ", var)
            next
         }
         
         if (gotLower + gotUpper == 2) {
            # We got some bounds
            if (!(lower < initial && initial < upper)) {
               warning("Check the bounds and ensure lower < initial < upper for ", var)
               next
            }
         }
         
         if (gotPrior) {
            if (prior <= 0) {
               warning("The prior standard deviation for ", var, " is not a positive value.")
               next
            }
         }
         
         # If the model has been estimated once then any initial state values named X and not X0
         # will have been changed to X0. Thus var must be changed as such.
         if (any(idx <- paste0(states, "0") %in% c(var, paste0(var, "0"))))
            var <- paste0(states[idx], "0")
         
         ParameterValues[var, "initial"] <<- initial
         ParameterValues[var, "lower"] <<- lower
         ParameterValues[var, "upper"] <<- upper
         ParameterValues[var, "prior"] <<- prior
         
         ParameterValues[var, "type"] <<- switch(gotInitial+gotLower+gotUpper+gotPrior,
                                                 "fixed",
                                                 "unused",
                                                 "ml",
                                                 "map")
      } # END : for (var in var.names) {
   } # END : setParameter(...) {
) # END : ctsm$methods()