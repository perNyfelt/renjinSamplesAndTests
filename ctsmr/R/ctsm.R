#' Continuous Time Stochastic Modelling for R
#' 
#' @author
#' Rune Juhl \email{ruju@@dtu.dk},
#' Henrik Madsen \email{hmad@@dtu.dk}
#' 
#' Maintainer: Rune Juhl \email{ruju@@dtu.dk}
#' @name ctsmr
#' @docType package
#' @aliases ctsmr package-ctsmr
#' @importFrom methods setRefClass
NULL

substituteCode <- function(e, pattern) {
   if (typeof(e) == "language") {
      ee <- as.list(e)
      for (i in seq_along(ee)) {
         if (typeof(ee[[i]]) %in% c("language","symbol")) {
            e[[i]] <- Recall(ee[[i]],pattern)
         }
      }
   } else if (as.character(e) == names(pattern)) {
      # We matched a leaf node with the pattern
      e <- pattern[[1]]
   }
   return(e)
}

checkDependence <- function(eqs, eqs_dep, eqno, notallowed=NULL) {
   notallowed <<- c(notallowed, eqno)
   eq <- eqs_dep[[eqno]]
   if (any(eq[notallowed])) { stop("Illegal cyclic dependence!") }
   if (any(eq)) {
      for (i in which(eq)) {
         eqs <- Recall(eqs, eqs_dep, i, notallowed)
         #eqs[[eqno]] <- substituteCode(eqs[[eqno]], eqs[i])
         eqs[[eqno]] <- codeWalker(eqs[[eqno]],insertEquations, eqs[i])
      }
   }
   return(eqs)
}

#' CTSM class
#'
#' CTSM Class (description)
#'
#' \code{ctsm} class is a reference class (details).
#'
#' @section Fields:
#' \itemize{
#'   \item N. Number of states
#'   \item M. Number of inputs
#'   \item S. Number of outputs
#'   \item NPARAM. Number of parameters
#' }
#' @section Methods:
#' \itemize{
#'   \item setOptions
#'   \item compile
#' }
#' @aliases ctsm
#' @export ctsm
#' @exportClass ctsm
ctsm <- setRefClass('ctsm',
      fields = list(
         N="numeric",M="numeric",S="numeric",NPARAM="numeric",Nest="numeric",
         alg.eqs="list",
         sys.eqs="list",
         obs.eqs="list",
         states="character",
         inputs="character",
         outputs="character",
         variance="list",
         linear="logical",
         uo="list",
         pars="character",
         libfile="character",
         ParameterValues="data.frame",
         nprior="numeric",
         PriorCorrelationMatrix="matrix",
         options="list"
      ),
      methods = list(
         initialize = function(...) {
            initFields(
               alg.eqs = list(), 
               sys.eqs = list(), 
               obs.eqs = list(),
               uo = list(), 
               options = .defaultoptions,
               linear = TRUE,
               nprior = 0,
               libfile = "",
               ParameterValues = data.frame(type    = numeric(0),
                                            initial = numeric(0),
                                            lower   = numeric(0),
                                            upper   = numeric(0),
                                            prior   = numeric(0))
            )
            #.self$addEquation # Force definition
            #.self$addEquations
         },
         finalize = function(...) {

         },
         #addEquation = addEquation2,
         #addEquations = addEquation2,
         setOptions = function(con=list(),default=FALSE) {
            'Use to change the default options.'
            
            if (default)
               options <<- .defaultoptions

            stopifnot(names(con) %in% names(.self$options))
            # solutionMethod         = "a",     
            # Method for solution propagation: subsampling, Adams, BDF
            if ("solutionMethod" %in% names(con)){
               con[["solutionMethod"]] <- switch(EXPR = con[["solutionMethod"]],
                         s           = 0,
                         sub         = 0,
                         subsampling = 0,
                         a           = 1,
                         adam        = 1,
                         adams       = 1,
                         nonstiff    = 1,
                         b           = 2,
                         bdf         = 2,
                         stiff       = 2)
            }
            options[names(con)] <<- con
         },
         
         estimate = function(data, firstorderinputinterpolation = FALSE, threads=1, debug=FALSE, verbose=FALSE) {
            'Estimate the model. Use threads >= 1 for parallel computation.'
            
            # Validate data
            data <- ValidateData(data, inputs, outputs, firstorderinputinterpolation)
            imat <- PackData(data, inputs, outputs)
            
            .self$AnalyseModel()
            
            parNames <- rownames(ParameterValues)
            # The initial states
            states0 <- sapply(states, paste0, "0")
            
            statedef <- ((states0 %in% parNames) + (states %in% parNames))
            
            if (any(statedef == 0))
               stop("Missing definition for initial state ", paste(states[statedef==0], collapse=","))
            
            if (!all(statedef == 1))
               stop("Multiple values for an initial state. Don't use both S and S0.")
            
            ii <- parNames %in% states
            if (any(ii)) {
               parNames[ii] <- states0[parNames[ii]]
               rownames(ParameterValues) <<- parNames
            }
            
            # Now all initial states should be called S0, V0 etc. pars should all be in the
            # ParameterValue
            if (!all(pars %in% parNames)) {
               stop("Missing initial values for ", paste(pars[!(pars %in% parNames)], collapse=", "), ".")
            }
            
            rho <- new.env(parent = baseenv())
            
            .self$options$sampletime <- as.numeric(!(attr(data, "ConstantSamplingTime")))
            .self$options$inputinterpolation <- attr(data, "FirstOrderInputInterpolation")
            
            assign(".nobs", as.integer(attr(imat, "Observations")), rho)
            assign(".nset", as.integer(length(data)), rho)
            
            ints <- genInts(ParameterValues,.self)
            assign(".ints", as.integer(ints), rho)
            assign(".nint", as.integer(length(ints)), rho)
            
            doublsold <- genDoubles(ParameterValues,.self)
            assign(".doubls", as.double(doublsold), rho)
            assign(".ndoub", as.integer(length(doublsold)), rho)
            
            assign(".tmat", as.double(attr(imat, "time")), rho)
            assign(".ntmat", as.integer(sum(attr(imat, "Observations"))), rho)
            
            assign(".imat", as.double(imat), rho)
            assign(".nimat", as.integer(length(imat)), rho)
            
            assign(".omat", as.double(0), rho)
            assign(".nomat", as.integer(1), rho)
                     
            assign(".trace", as.double(numeric((3*Nest+2)*options$maxNumberOfEval)), rho)
            assign(".info", as.integer(0), rho)
            
            assign(".cpus", as.integer(0),rho)
            assign(".threads", as.integer(threads),rho)
            
            if (debug) {
               rawinput = list(
                  tmat = as.double(attr(imat, "time")),
                  imat = as.double(imat),
                  ints = as.integer(ints),
                  doubls = as.double(doublsold))
            }
            
            CallFunction("ctsmdriver", .self, rho, verbose = verbose)
            
            raw <- mget(c(".ints",".doubls",".trace",".threads",".cpus",".info"), envir=rho)
            
            # Prepare a return object
            z <- list()
            
            if (raw$.info != 0) {
               warning(sprintf("Did not converge. %s", infocode(raw$.info)))
            } else {
               
#                pars2 <- as.logical(c(estState$est,estPar$est))
#                names(pars2) <-c(estState$name,estPar$name)
#                pars2 <- pars2[match(pars,names(pars2))]
#                
#                I <- seq(1,sum(pars2))
#                J <- which(pars2)
               
               # Which parameters where estimated?
               pars2 <- ParameterValues[pars,"type"] != "fixed"
               names(pars2) <- pars
               
               I <- seq(1,Nest)
               J <- which(pars2)
               
               parnames <- pars[pars2]
               f        = raw$.doubls[4]
               fpen     = raw$.doubls[(NPARAM+5)*NPARAM+5]
               fprior   = raw$.doubls[(NPARAM+5)*NPARAM+6]
               
               z <- list(         
                  itr      = raw$.ints[2],
                  neval    = raw$.ints[3],
                  
                  detH     = raw$.doubls[3],
                  f        = f,
                  fpen     = fpen,
                  fprior   = fprior,
                  loglik   = fpen - f,
                  
                  dL       = raw$.doubls[4+I],
                  dpen     = raw$.doubls[NPARAM+4+I],
                  
                  corr     = addNames(matrix(
                     raw$.doubls[seq(2*NPARAM+5,(NPARAM+2)*NPARAM+4)],
                     NPARAM, NPARAM, byrow=F)[I,I],
                                      parnames),
                  
                  pt       = addNames(raw$.doubls[(NPARAM+2)*NPARAM+4+I], parnames),
                  sd       = addNames(raw$.doubls[(NPARAM+3)*NPARAM+4+I], parnames),
                  t        = addNames(raw$.doubls[(NPARAM+4)*NPARAM+4+I], parnames),
                  
                  xm       = addNames(raw$.doubls[(2*NPARAM+6)*NPARAM+12+seq(1,NPARAM)], pars),
                  xmax     = addNames(raw$.doubls[(2*NPARAM+7)*NPARAM+12+I], parnames),
                  xmin     = addNames(raw$.doubls[(2*NPARAM+8)*NPARAM+12+I], parnames),
                  
                  trace    = matrix(raw$.trace, ncol=3*Nest+2, byrow=TRUE),
                  
                  estimated = pars2
                  )
            }
            
            z$threads  <- raw$.threads
            z$cpus     <- raw$.cpus
            z$info     <- raw$.info
            z$message  <- infocode(raw$.info)
            
            if (debug) {
               z$rawinput <- rawinput
               z$env <- rho
            }
            
            z$model <- .self$copy(shallow=TRUE)
            
            z$data <- data
            
            class(z) <- "ctsmr"
            rm(rho)
            return(z)
         },
         
         show = function() {
            'Prints info about the model.'
            
            n <- length(states)
            m <- length(inputs)
            s <- length(outputs)
            
            if (n+m+s == 0) {
               cat("Empty CTSM-R model.\n")  
            } else {
               
               cat(paste0(ifelse(linear, "Linear", "Non-linear"), ' state space model with ', n, ' state', ifelse(n>1,"s",""), ", ",
                          s, ' output', ifelse(s>1, "s", ""), " and ", m, " input", ifelse(m>1,"s",""),"\n"))
               
               if (n) {
                  cat("\nSystem equations:\n")
                  cat(paste0("\t", 
                             paste(lapply(states, function(x) paste(deparse(sys.eqs[[x]]$input), collapse="\n\t")), collapse="\n\t"),"\n"))
               }
               
               if (s) {
                  cat("\nObservation equations:\n")
                  cat(paste0("\t", 
                             paste(lapply(outputs, function(x) paste(deparse(obs.eqs[[x]]$input), collapse="\n\t")), collapse="\n\t"),"\n"))
               }
               
               if (m) {
                  cat("\nInputs: ", paste(inputs, collapse=", "), "\n")
               } else
                  cat("\nNo inputs.\n")
               
               if (length(pars)) {
                  cat("\nParameters: ", paste(pars, collapse=", "), "\n")
               }
            }
         }
      )
)

ctsm$methods(
   addSystem = function(form) {
      'Adds one or more stochastic differential equations to the model. The diffusion processes must be called dw{1,2,...}'
      res <- CheckSystemEquation(form, .self)
      CheckExistance(names(res), .self)
      sys.eqs <<- c(sys.eqs, res)
      states <<- c(states, names(res))
   }
)

ctsm$methods(
   addObs = function(form) {
      'Adds one or more observation equations to the model.'
      res <- CheckObservationEquation(form, .self)
      CheckExistance(names(res), .self)
      obs.eqs <<- c(obs.eqs, res)
      outputs <<- c(outputs, names(res))
   }
)

ctsm$methods(
   addInput = function(...) {
      'Specifies which variables are inputs.'
      
      # Intercept the call
      cl <- match.call()
      
      # Any arguments at all?
      if (length(cl) == 1) { stop("Need at least one parameter") }
      
      input.vars <- as.character(as.list(cl[-1]))
      valid <- grepl(pattern="^[[:alpha:]][[:alnum:]]*$",x=input.vars)
      if (!all(valid)) {
         if (sum(!valid) == 1)
            text <- "Illegal input variable name:"
         else
            text <- "Illegal input variable names:"
         stop( paste(text, paste(input.vars[!valid], collapse=", ")) )
      }
      
      # Ok! All input strings are valid so far
      
      # Check if the inputs exist as states
      chk.states <- is.na(match(states, input.vars))
      if (!all(chk.states)) {
         if (sum(!chk.states) == 1)
            text <- "Input variable already exist as a state:"
         else
            text <- "Input variables already exist as states:"
         stop( paste(text, paste(states[!chk.states], collapse=", ")) )
      }
      
      inputs <<- unique(c(inputs, input.vars))
   }
)

ctsm$methods(
   setVariance = function(form) {
      'Sets the variance or covariance of the outputs.'
      
      if (class(form) != "formula") { stop("Argument must be a formula.") }
      
      res <- ParseOutputVariance(form , outputs)
      
      if (is.null(res)) { stop("Failed to understand the provided variance term.") }
      
      variance[[length(variance)+1]] <<- res
   }
)

ctsm$methods(
   AnalyseModel = function() {
      'Analyses the structure of the model.'
      
      if (length(states) == 0) { stop("Need at least one system equation") }
      if (length(outputs) == 0) { stop("Need at least one observation equation") }
      
      for (state in states)
         sys.eqs[[state]]$processed <<- sys.eqs[[state]]$input
      
      for (output in outputs)
         obs.eqs[[output]]$processed <<- obs.eqs[[output]]$input
      
      # First expand the algebraic equations
      if (length(.self$alg.eqs) != 0) {
         alg.eqs.processed <- ExpandEquations(.self$alg.eqs)$parsed
         # Insert the equations 
         for (i in states)
            sys.eqs[[i]]$processed <<- lapply(sys.eqs[[i]]$input, 
                                              function(x) { codeWalker(x, insertEquations, alg.eqs.processed) })
         for (i in outputs)
            obs.eqs[[i]]$processed <<- lapply(obs.eqs[[i]]$input, 
                                              function(x) { codeWalker(x, insertEquations, alg.eqs.processed) })
      }
      
      for (i in states) {
         res <- CheckSystemEquation(sys.eqs[[i]]$processed, .self)[[1]]
         
         sys.eqs[[i]]$parsed <<- res$parsed
         sys.eqs[[i]]$invariant <<- res$invariant
         
         res <- IsLinear(sys.eqs[[i]]$parsed$dt, .self)
         
         # Kinda nasty hack to handle cases which are linear in states and inputs, but with inputs or "t" in the diffusion term
         res$linear <- res$linear && (!any(c(inputs,"t") %in% getPars(sys.eqs[[i]]$parsed[names(sys.eqs[[i]]$parsed) != "dt"])))
         
         sys.eqs[[i]]$linear <<- res$linear
         sys.eqs[[i]]$f <<- sys.eqs[[i]]$parsed$dt
         
         if (res$linear) {
            # L case
            sys.eqs[[i]]$A <<- res$A
            sys.eqs[[i]]$B <<- res$B
         } else
            linear <<- FALSE      
      }
      
      for (i in outputs) {
         res <- CheckObservationEquation(obs.eqs[[i]]$processed, .self)[[1]]
         
         obs.eqs[[i]]$parsed <<- res$parsed
         obs.eqs[[i]]$invariant <<- res$invariant
         
         res <- IsLinear(obs.eqs[[i]]$parsed, .self)
         obs.eqs[[i]]$linear <<- res$linear
         obs.eqs[[i]]$h <<- obs.eqs[[i]]$parsed
         
         if (res$linear) {
            # L case
            obs.eqs[[i]]$C <<- res$A
            obs.eqs[[i]]$D <<- res$B
         } else
            linear <<- FALSE      
      }
      
      # Check for inputs or "t" in the variance symbols
      variance_symbols <- lapply(variance, function(x) getPars(x$var))
      # If any inputs or "t" in variance_symbols then the model is nonlinear (LTV in fact, but hey..)
      linear <<- linear && !any(c(inputs, "t") %in% variance_symbols)
      
      if (linear) {
         options$solutionMethod <<- 0
         options$numberOfSubsamples <<- 1
         options$nIEKF <<- 1
         # Set the type = 0 - ONLY FOR LTI!
         options$type <<- 0
      } else
         options$type <<- 2
      
      N <<- length(states)
      M <<- length(inputs)
      S <<- length(outputs)
      
      # Find all parameters
      
      pars <<- unique(unlist(c(lapply(sys.eqs[states], function (x) { getPars(x$parsed) }),
                               lapply(obs.eqs[outputs], function (x) { getPars(x$parsed) }),
                               variance_symbols )))
      
      pars <<- unique(c(gsub("(.+)","\\10",states), 
                       sort(setdiff(pars,c(states,inputs)))))
      
      if ("t" %in% pars) {
         pars <<- pars[pars != "t"]
      }
      
      # Is the initial covariance matrix k*I? If so add scaleVariance as a parameter.
      if (options$InitialVarianceScaledIdentity) {
         pars <<- c(pars, "scaleVariance")
         if (!("scaleVariance" %in% rownames(ParameterValues))) {
            .self$setParameter(scaleVariance = c(log(1e2),-18,log(1e4)))
         }
      }
      
      NPARAM <<- length(pars)
   }
)


#calltoFunction <- function(n,m,fcns,env,xpvar=NULL,xpval=NULL,xmvar=NULL,xmval=NULL,uovar=NULL,uoval=NULL) {

calltoFunction <- function(fcns,env,xpvar=NULL,xmvar=NULL,uovar=NULL) {
   asg <- as.name("<-")
   brc <- as.name("[")
   mat <- as.name("mat")

   f <- function(xp,xm,uo,t) { NULL }
   i = 1
   if (!(is.null(xpvar)))
      for (j in 1:length(xpvar))
         body(f)[[i<-i+1]] <- as.call(list(asg,as.name(xpvar[j]), as.call(list(brc,as.name("xp"),as.numeric(j)))))

   if (!(is.null(xmvar)))
      for (j in 4:length(xmvar))
         body(f)[[i<-i+1]] <- as.call(list(asg,as.name(xmvar[j]), as.call(list(brc,as.name("xm"),as.numeric(j)))))

   if (!(is.null(uovar)))
      for (j in 1:length(uovar))
         body(f)[[i<-i+1]] <- as.call(list(asg,as.name(uovar[j]), as.call(list(brc,as.name("uo"),as.numeric(j)))))
   
   tmp <- as.call(list(as.name("array"),0.,length(fcns)))
   body(f)[[i<-i+1]] <- as.call(list(asg,mat,tmp))

   for (j in 1:length(fcns)) {
      if (!(is.null(fcns[[j]]) || fcns[[j]] == 0)) {
         tmp <- as.call(list(brc,mat,as.numeric(j)))
         body(f)[[i<-i+1]] <- as.call(list(asg,tmp,fcns[[j]]))
      }
   }
   body(f)[[i+1]] <- quote(mat)

   environment(f) <- env
   return(f)
}

ctsm$methods(
   pred = function(res, newdata, opt, k=NULL, covariance=FALSE) {
      'General prediction, smoothing, filtering method.'
      
      if (covariance == TRUE && opt == "p") {
         warning("The full covariance cannot be obtained when predicting.")
         covariance <- FALSE
      }
      
      rho = new.env(parent=baseenv())
      
      has.time <- newdata$data[["t"]]
      
      Nt <- length(has.time)
      
      nparam <- NPARAM
      
      .self$options$sampletime <- newdata$sampletime
      .self$options$inputinterpolation <- newdata$interpolation
      
      assign(".nobs", as.integer(Nt), rho)
      assign(".nset", 1L, rho)
      
      ints <- genInts(NULL, .self, pred=TRUE)
      assign(".ints", as.integer(ints), rho)
      assign(".nint", as.integer(length(ints)), rho)
      
#       stateval <- list(name = pars[1:N],initval = res$xm[pars[1:N]])
#       parval   <- list(name = pars[-c(1:N)],initval = res$xm[pars[-c(1:N)]])
      
#       doubls <- genDoubles(stateval,parval,.self,pred=TRUE)
      doubls <- genDoubles(data.frame(initial=res$xm), .self, pred=TRUE)
      
      assign(".doubls", as.double(doubls), rho)
      assign(".ndoub", as.integer(length(doubls)), rho)
      
      assign(".tmat", as.double(has.time), rho)
      assign(".ntmat", as.integer(length(has.time)), rho)
            
      imat <- newdata$data[c(inputs,outputs)]
      imat <- c(t(imat))
      
      assign(".imat", as.double(imat), rho)
      assign(".nimat", as.integer(length(imat)), rho)
      
      if (opt == "f") {
         Ncomat <- 2*N + as.double(covariance)*N*(N-1)/2
         omat <- numeric(Nt * Ncomat)
         omat[1] <- as.double(covariance)
         assign(".job", as.integer(-2), rho)
      } else if (opt == "s") {
         Ncomat <- 2*N + as.double(covariance)*N*(N-1)/2 +
                   2*S + as.double(covariance)*S*(S-1)/2
         omat <- numeric(Nt * Ncomat)
         
         omat[1] <- as.double(covariance*3)
         assign(".job", as.integer(-1), rho)
      } else if (opt == "p") {
         Ncomat <- 2 * N + 2 * S
         omat <- numeric(Nt * Ncomat)
         assign(".job", as.integer(k), rho)
      }

      assign(".omat", as.double(omat), rho)
      assign(".nomat", as.integer(length(omat)), rho)
      # TODO: Wrong size. Should be #estimated vars
#       Nest <<- sum(estState$est)+sum(estPar$est)
      
      assign(".info", as.integer(0), rho)
      
      fun <- getNativeSymbolInfo("ctsmllike",lib)
      ans <- .Call(fun,rho)
      
      info <- get(".info", rho)
      
      if (info != 0)
         stop(infocode(info))
      
      #omat <- matrix(get(".omat", rho),ncol=2*.self$N+cov*.self$N*(.self$N-1)/2)
      omat <- get(".omat", rho)
      rm(".omat", envir=rho)
      
      omat <- matrix(omat, ncol=Ncomat, byrow=TRUE)
      
      if (opt == "f") { # Filtering
         
         if (covariance) {
            covstate <- ToCovariance(N, omat[,-c(1:N)])
            
            return(list(
               state=omat[,1:N],
               covariance=covstate)
                   )
         } else {
            # Otherwise return the standard deviations
            return(list(
               state=omat[,1:N],
               sd=omat[,-c(1:N)]))
         }
      } else if (opt == "s") { # Pure simulation
         
         state.value <- omat[, c(1:(k<-N))]
         state.var   <- omat[, c((k+1):(k<-k+N+as.double(covariance)*N*(N-1)/2))]
         
         obs.value   <- omat[, c((k+1):(k<-k+S))]
         obs.var     <- omat[, c((k+1):(k<-k+S+as.double(covariance)*S*(S-1)/2))]
         
         if (covariance) {
            state.var <- ToCovariance(N, state.var)
            obs.var   <- ToCovariance(S, obs.var)
            
            return(list(
               state = list(value=state.value, covariance=state.var),
               obs   = list(value=obs.value, covariance=obs.var))
                   )
         } else {
            # Otherwise return the standard deviations
            return(list(
               state = list(value=state.value, sd=state.var),
               obs   = list(value=obs.value, sd=obs.var))
                   )
         }
         
      } else if (opt == "p") { # k-step prediction
         return(list(
            state  = list(pred = omat[,1:.self$N],
                                   sd   = omat[,(.self$N+1):(2*.self$N)]),
            output = list(
               pred = omat[,(2*.self$N+1):(2*.self$N+.self$S)],
               sd   = omat[,(2*.self$N+.self$S+1):(2*.self$N+2*.self$S)])
            )
                )
      }
   }
   )

ctsm$methods(
   conditionallikelihood = function(fit, which, maxpts = 100, data) {
      'Calculates a conditional likelihood.'
      
      stopifnot(is.character(which) || is.name(which))
      
      if (!(which %in% pars)) {
         stop(which, " is not a paramater in the model.")
      }
      
      if (ParameterValues[which, "type"] == "fixed") {
         stop(which, " is a fixed parameter.")
      }
      
      rho = new.env(parent=baseenv())

      datanames <- names(data$data)
      
      timename <- c("t","time")
      stopifnot(any(has.time <- c("t","time") %in% datanames))

      has.time <- data$data[[timename[has.time]]]
      
      
      nparam <- NPARAM
      
      .self$options$sampletime <- data$sampletime
      .self$options$inputinterpolation <- data$interpolation
      
      assign(".nobs", as.integer(length(has.time)), rho)
      assign(".nset", as.integer(1), rho)
      
#       ints <- genInts(estState,estPar,.self)
#       
#       ints[9+(1:nparam)] <- 0
      
      newPar <- ParameterValues
      newPar[names(fit$xm), "initial"] <- fit$xm
      
      ints <- genInts(newPar,.self)
      assign(".ints", as.integer(ints), rho)
      assign(".nint", as.integer(length(ints)), rho)
      
      doubls <- genDoubles(newPar,.self)
      
      
      # Transform all estimated parameters
      # log( (x - xmin) / (xmax - x))
      newPar <- newPar[pars, ]
      newPar <- newPar[newPar[,"type"] == "ml", ]
      
      xi <- which(rownames(newPar) == which)
      
      xx <- newPar$initial
      xxmin <- newPar$lower
      xxmax <- newPar$upper
      
      doubls[(2*nparam+7)*nparam+12+(1:nparam)] <- c(xxmax,numeric(nparam-length(xx)))
      doubls[(2*nparam+8)*nparam+12+(1:nparam)] <- c(xxmin,numeric(nparam-length(xx)))
      
      assign(".doubls", as.double(doubls), rho)
      assign(".ndoub", as.integer(length(doubls)), rho)
      
      assign(".tmat", as.double(has.time), rho)
      assign(".ntmat", as.integer(length(has.time)), rho)
      
      imat <- data$data[c(inputs,outputs)]
      imat <- c(t(imat))
      
      assign(".imat", as.double(imat), rho)
      assign(".nimat", as.integer(length(imat)), rho)
      
      assign(".job", as.integer(0), rho)
      
      assign(".omat", as.double(0), rho)
      assign(".nomat", as.integer(1), rho)
      
      assign(".nmisst", as.integer(0), rho)
      assign(".f", as.double(0), rho)
      assign(".fpen", as.double( 0 ), rho)
      assign(".fprior", as.double(0), rho)
      
      assign(".info", as.integer(0), rho)
      
      fun <- getNativeSymbolInfo("ctsmllike2", lib)
      
      f <- numeric(maxpts)
      fpen <- numeric(maxpts)
      
      l <- ParameterValues[which, "lower"]
      u <- ParameterValues[which, "upper"]
      d <- (u-l) * 0.05
      
      xseq <- seq(l+d, u-d, length.out=maxpts)
      
      assign(".nx", as.integer(length(xx)), rho)
      
      for (i in 1:maxpts) {
         
         xx[xi] <- xseq[i]
         
         xt <- log( (xx - xxmin) / (xxmax - xx))
         
         assign(".x", as.double(xt), rho)
         
         
         ans <- .Call(fun,rho)
         
         info <- get(".info", rho)
         
         if (info != 0)
            stop(infocode(info))
         
         f[i] <- get(".f", rho)
         fpen[i] <- get(".fpen", rho)
      
      }
      
      return(list(f=f,fpen=fpen,xseq=xseq))
      
   }
)

ctsm$methods(
   profilelikelihood = function(fit, which, data, maxpts = 100, nsd = 5, threads=1) {
      'Calculates a profile likelihood.'
      
      stopifnot(is.character(which) || is.name(which))
      
      if (!(which %in% pars)) {
         stop(which, " is not a paramater in the model.")
      }
      
      if (ParameterValues[which, "type"] == "fixed") {
         stop(which, " is a fixed parameter.")
      }
      
      rho = new.env(parent=baseenv())
      
      datanames <- names(data$data)
      
      timename <- c("t","time")
      stopifnot(any(has.time <- c("t","time") %in% datanames))
      
      has.time <- data$data[[timename[has.time]]]
      
      
      nparam <- NPARAM
      
      .self$options$sampletime <- data$sampletime
      .self$options$inputinterpolation <- data$interpolation
      
      
      assign(".nobs", as.integer(length(has.time)), rho)
      assign(".nset", as.integer(1), rho)
         
      newPar <- ParameterValues
      newPar[names(fit$xm), "initial"] <- fit$xm
      newPar[which, "type"] <- "fixed"
      
      
      ints <- genInts(newPar,.self)
      assign(".ints", as.integer(ints), rho)
      assign(".nint", as.integer(length(ints)), rho)
      
      doubls <- genDoubles(newPar,.self)
         
      assign(".doubls", as.double(doubls), rho)
      assign(".ndoub", as.integer(length(doubls)), rho)
      
      assign(".tmat", as.double(has.time), rho)
      assign(".ntmat", as.integer(length(has.time)), rho)
      
      imat <- data$data[c(inputs,outputs)]
      imat <- c(t(imat))
      
      assign(".imat", as.double(imat), rho)
      assign(".nimat", as.integer(length(imat)), rho)
      
      assign(".job", as.integer(0), rho)
      
      assign(".omat", as.double(0), rho)
      assign(".nomat", as.integer(1), rho)
      
      assign(".trace", as.double(numeric((3*Nest+2)*options$maxNumberOfEval)), rho)
      assign(".info", as.integer(0), rho)
      
      assign(".cpus", as.integer(0),rho)
      assign(".threads", as.integer(threads),rho)
      
      fun <- getNativeSymbolInfo("ctsmdriver", lib)
      
      f <- fpen <- itr <- neval <- numeric(maxpts)
      xmm <- matrix(nrow=maxpts, ncol=NPARAM)
      colnames(xmm) <- pars
      
      l <- fit$xm[which] - nsd*fit$sd[which]
      u <- fit$xm[which] + nsd*fit$sd[which]
      
      xseq <- seq(l, u, length.out=maxpts)
      
      #xseq <- c(seq(0,res$xm[which],length.out=100)[-100], seq(res$xm[which], 0.02, length.out=50))
      
      istart <- which.min((xseq - newPar[which, "initial"])^2)
      
      cat("Profiling started")
      
      for (i in istart:maxpts) {
         
         newPar[which, "initial"] <- xseq[i]
         
         doubls <- genDoubles(newPar,.self)
         
         assign(".doubls", as.double(doubls), rho)
         assign(".ints", as.integer(ints), rho)
         
         cat(paste("i =", i))
         
         ans <- .Call(fun,rho)
         
         info <- get(".info", rho)
         
         if (info != 0) {
            warning(infocode(info))
            next
         }
                  
         tmp <- get(".doubls", rho)
         xmm[i,] <- tmp[(2*NPARAM+6)*NPARAM+12+seq(1,NPARAM)]
         newPar[pars, "initial"] <- xmm[i,]
         
         f[i]        <- tmp[4]
         fpen[i]     <- tmp[(NPARAM+5)*NPARAM+5]
         #fprior[i]   = tmp[(NPARAM+5)*NPARAM+6]
         
         tmp <- get(".ints", rho)
         
         itr[i] <- tmp[2]
         neval[i] <- tmp[3]
         
      }
      
      newPar[pars, "initial"] <- xmm[istart,]
      
      for (i in (istart-1):1) {
                  
         newPar[which, "initial"] <- xseq[i]
         
         doubls <- genDoubles(newPar,.self)
         
         assign(".doubls", as.double(doubls), rho)
         assign(".ints", as.integer(ints), rho)
         
         cat(paste("i =", i))
         
         ans <- .Call(fun,rho)
         
         info <- get(".info", rho)
         
         if (info != 0) {
            warning(infocode(info))
            next
         }
         
         tmp <- get(".doubls", rho)
         xmm[i,] <- tmp[(2*NPARAM+6)*NPARAM+12+seq(1,NPARAM)]
         newPar[pars, "initial"] <- xmm[i,]
         
         f[i]        <- tmp[4]
         fpen[i]     <- tmp[(NPARAM+5)*NPARAM+5]
         #fprior[i]   = tmp[(NPARAM+5)*NPARAM+6]
         
         tmp <- get(".ints", rho)
         
         itr[i] <- tmp[2]
         neval[i] <- tmp[3]
         
      }
      
      return(list(f=f,fpen=fpen,xseq=xseq,itr=itr,neval=neval))
      
   }
)

#setMethod("$<-", "ctsm", 
#          function(x, name, value) {
#             
#             if (name %in% c("sys.eqs","obs.eqs"))
#                stop(sprintf("%s cannot be altered.", name), domain=NA)
#             
#             x[[name]] <- value
#             
#             return(x)
#          }, sealed=TRUE)
#