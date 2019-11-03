#' Predict method for CTSM fits
#' @param object Object of class 'ctsmr'
#' @param n.ahead The number of steps ahead for which prediction is required.
#' @param covariance Should the full covariance matrix be outputted.
#' @param newdata An optional data frame with new data.
#' @param firstorderinputinterpolation If TRUE (default FALSE) first order linear interpolation of the inputs between samples.
#' @param x0 initial state
#' @param vx0 initial state covariance
#' @param ... Unused.
#' @export
predict.ctsmr <- function(object, n.ahead = 1, covariance = FALSE, newdata = NULL, firstorderinputinterpolation = FALSE, x0 = NULL, vx0 = NULL, ...) {
   
   ValidateResultObject(object)
   
   if (missing(newdata) || is.null(newdata)) {
      # Use data used while fitting
      newdata <- object$data
      firstorderinputinterpolation <- attr(object$data, "FirstOrderInputInterpolation")
   }
   
   m <- object$model
   
   outputsize <- 2L*(m$N + m$S) + ifelse(covariance, m$N*(m$N-1L)/2L + m$S*(m$S-1L)/2L, 0L)
   
   e <- PrepareEnv(object, 
                   newdata,
                   firstorderinputinterpolation = firstorderinputinterpolation,
                   job = n.ahead,
                   outputsize = outputsize,
                   covariance = ifelse(covariance, 3L, 0L),
                   x0 = x0,
                   vx0 = vx0
   )
   
   CallFunction("ctsmllike", m, e)
   
   if ((info <- get(".info", e)) != 0)
      stop(infocode(info))
   
   omat <- matrix(get(".omat", e), ncol=outputsize, byrow=TRUE)
   
   to_named_dataframe <- function(mat, vnames, covariance = FALSE) {
      if (!covariance) {
         mat <- as.data.frame(mat)
         colnames(mat) <- vnames
      } else {
         mat <- ToCovariance(length(vnames), mat, vnames)
      }
      mat
   }
   
   state.pred <- to_named_dataframe(omat[, 1L:m$N ], m$states)
   
   state.var  <- to_named_dataframe(omat[, m$N + seq.int(1L, cpos <- m$N * ifelse(covariance, (m$N + 1L)/2L, 1L)) ],
                                    m$states, covariance)
   cpos <- cpos + m$N
   
   output.pred <- to_named_dataframe(omat[, cpos + seq.int(1L, m$S) ], m$outputs)
   cpos <- cpos + m$S
   
   output.var  <- to_named_dataframe(omat[, cpos + seq.int(1L, m$S * ifelse(covariance, (m$S + 1L)/2L, 1L)) ],
                                     m$outputs, covariance)
   
   
   state = list(state.pred, state.var)
   names(state) <- c("pred", ifelse(covariance, "var", "sd"))
   
   output = list(output.pred, output.var)
   names(output) <- c("pred", ifelse(covariance, "var", "sd"))
   
   if (is.data.frame(newdata)) {
      pred <- list(state=state, output=output)      
   } else {
      pred <- vector(mode="list", length=length(newdata))
      Nprev <- 0L
      for (i in seq_along(newdata)) {
         idx <- seq.int(1L, nrow(newdata[[i]])) + Nprev
         Nprev <- Nprev + nrow(newdata[[i]])
         
         pred[[i]] <- list(state  = list(pred = state.pred[idx, , drop=FALSE],
                                         if (covariance) {
                                            state.var[, , idx]
                                         } else {
                                            state.var[idx, , drop=FALSE]
                                         }
         ),
         output = list(pred = output.pred[idx, , drop=FALSE],
                       if (covariance) {
                          output.var[, , idx]
                       } else {
                          output.var[idx, , drop=FALSE]
                       }
         )
         )
         names(pred[[i]]$state) <- c("pred",ifelse(covariance,"var","sd"))
         names(pred[[i]]$output) <- c("pred",ifelse(covariance,"var","sd"))
      }
      names(pred) <- names(newdata)
   }
   
   #    colnames(pred$state$pred) <- m$states
   #    colnames(pred$state$sd) <- m$states
   #    
   #    colnames(pred$output$pred) <- m$outputs
   #    colnames(pred$output$sd) <- m$outputs
   
   pred
}

#' Optimal filter method for CTSM fits
#' @param object Object of class 'ctsmr'
#' @param newdata An optional data frame with new data.
#' @param firstorderinputinterpolation If TRUE (default FALSE) first order linear interpolation of the inputs between samples.
#' @param covariance If TRUE (default FALSE) the full covariance matrix will be returned.
#' @param x0 initial state
#' @param vx0 initial state covariance
#' @param ... Unused.
#' @export
filter.ctsmr <- function(object, newdata = NULL, firstorderinputinterpolation = FALSE, covariance = FALSE, x0 = NULL, vx0 = NULL, ...) {
   
   ValidateResultObject(object)
   
   if (missing(newdata) || is.null(newdata)) {
      # Use data used while fitting
      newdata <- object$data
      firstorderinputinterpolation <- attr(object$data, "FirstOrderInputInterpolation")
   }
   
   m <- object$model
   
   outputsize <- 2L*(m$N) + ifelse(covariance, m$N*(m$N-1L)/2L, 0L)
   
   e <- PrepareEnv(object, 
                   newdata,
                   firstorderinputinterpolation = firstorderinputinterpolation,
                   job = -2L,
                   outputsize = outputsize,
                   covariance = ifelse(covariance, 3L, 0L),
                   x0 = x0,
                   vx0 = vx0
   )
   
   CallFunction("ctsmllike", m, e)
   
   if ((info <- get(".info", e)) != 0)
      stop(infocode(info))
   
   omat <- matrix(get(".omat", e), ncol=outputsize, byrow=TRUE)
   
   to_named_dataframe <- function(mat, vnames, covariance = FALSE) {
      if (!covariance) {
         mat <- as.data.frame(mat)
         colnames(mat) <- vnames
      } else {
         mat <- ToCovariance(length(vnames), mat, vnames)
      }
      mat
   }
   
   state.pred <- to_named_dataframe(omat[, 1L:m$N ], m$states)
   
   state.var  <- to_named_dataframe(omat[, m$N + seq.int(1L, cpos <- m$N * ifelse(covariance, (m$N + 1L)/2L, 1L)) ],
                                    m$states, covariance)
   
   state = list(state.pred, state.var)
   names(state) <- c("filtered", ifelse(covariance, "var", "sd"))
   
   if (is.data.frame(newdata)) {
      pred <- state      
   } else {
      pred <- vector(mode="list", length=length(newdata))
      Nprev <- 0L
      for (i in seq_along(newdata)) {
         idx <- seq.int(1L, nrow(newdata[[i]])) + Nprev
         Nprev <- Nprev + nrow(newdata[[i]])
         
         pred[[i]] <- list(filtered = state.pred[idx, , drop=FALSE],
                           if (covariance) {
                              state.var[, , idx]
                           } else {
                              state.var[idx, , drop=FALSE]
                           }
         )
         names(pred[[i]]) <- c("filtered", ifelse(covariance,"var","sd"))
      }
      names(pred) <- names(newdata)
   }
   
   #    colnames(pred$state$pred) <- m$states
   #    colnames(pred$state$sd) <- m$states
   #    
   #    colnames(pred$output$pred) <- m$outputs
   #    colnames(pred$output$sd) <- m$outputs
   
   pred
   
}

#' Smoothing filter method for CTSM fits
#' @param object Object of class 'ctsmr'
#' @param newdata An optional data frame with new data.
#' @param firstorderinputinterpolation  If TRUE (default FALSE) first order linear interpolation of the inputs between samples.
#' @param ... Unused.
#' @export
smooth.ctsmr <- function(object, newdata = NULL, firstorderinputinterpolation = FALSE, ...) {
   
   ValidateResultObject(object)
   
   if (missing(newdata) || is.null(newdata)) {
      # Use data used while fitting
      newdata <- object$data
      firstorderinputinterpolation <- attr(object$data, "FirstOrderInputInterpolation")
   }
   
   m <- object$model
   e <- PrepareEnv(object=object, data=newdata,
                   firstorderinputinterpolation = firstorderinputinterpolation, 
                   job = -3L,
                   outputsize = 2L*(m$N+m$N*(m$N+1L)/2L))
   
   CallFunction("ctsmllike", object$model, e)
   
   info <- get(".info", e)
   
   if (info != 0L)
      stop(infocode(info))
   
   omat <- get(".omat", e)
   rm(".omat", envir=e)
   
   omat <- matrix(omat, ncol=2L*m$N, byrow=TRUE)
   
   if (is.data.frame(newdata)) {
      ans <- list(
         state = omat[1:nrow(newdata$data),    1:m$N ],
         sd    = omat[1:nrow(newdata$data), -c(1:m$N)])      
   } else {
      ans <- vector(mode="list", length=length(newdata))
      Nprev <- 0
      for (i in seq_along(newdata)) {
         idx <- seq(1, nrow(newdata[[i]])) + Nprev
         Nprev <- Nprev + nrow(newdata[[i]])
         ans[[i]] <- list(state = omat[idx, 1:m$N],
                          sd    = omat[idx, (m$N+1):(2*m$N)])
      }
      names(ans) <- names(newdata)
   }
   
   return( ans )
}

#' Simulate method for CTSM fits
#' @param object Object of class 'ctsmr'
#' @param nsim number of realisations. Defaults to 1.
#' @param seed seed
#' @param newdata An optional data frame with new data.
#' @param firstorderinputinterpolation If TRUE (default FALSE) first order linear interpolation of the inputs between samples.
#' @param covariance If TRUE (default FALSE) the full covariance matrix will be returned.
#' @param x0 initial state
#' @param vx0 initial state covariance
#' @param ... Unused.
#' @export
#' @importFrom stats simulate
simulate.ctsmr <- function(object, nsim = 1L, seed = NULL, newdata = NULL, firstorderinputinterpolation = FALSE, covariance = FALSE, x0 = NULL, vx0 = NULL, ...) {

   ValidateResultObject(object)
   
   if (missing(newdata) || is.null(newdata)) {
      # Use data used while fitting
      newdata <- object$data
      firstorderinputinterpolation <- attr(object$data, "FirstOrderInputInterpolation")
   }
   
   m <- object$model
   
   outputsize <- 2L*(m$N + m$S) + ifelse(covariance, m$N*(m$N-1L)/2L + m$S*(m$S-1L)/2L, 0L)
   
   e <- PrepareEnv(object, 
                   newdata,
                   firstorderinputinterpolation = firstorderinputinterpolation,
                   job = -1L,
                   outputsize = outputsize,
                   covariance = ifelse(covariance, 3L, 0L),
                   x0 = x0,
                   vx0 = vx0
   )
   
   CallFunction("ctsmllike", m, e)
   
   if ((info <- get(".info", e)) != 0)
      stop(infocode(info))
   
   omat <- matrix(get(".omat", e), ncol=outputsize, byrow=TRUE)
   
   to_named_dataframe <- function(mat, vnames, covariance = FALSE) {
      if (!covariance) {
         mat <- as.data.frame(mat)
         colnames(mat) <- vnames
      } else {
         mat <- ToCovariance(length(vnames), mat, vnames)
      }
      mat
   }
   
   state.pred <- to_named_dataframe(omat[, 1L:m$N ], m$states)
   
   state.var  <- to_named_dataframe(omat[, m$N + seq.int(1L, cpos <- m$N * ifelse(covariance, (m$N + 1L)/2L, 1L)) ],
                                    m$states, covariance)
   cpos <- cpos + m$N
   
   output.pred <- to_named_dataframe(omat[, cpos + seq.int(1L, m$S) ], m$outputs)
   cpos <- cpos + m$S
   
   output.var  <- to_named_dataframe(omat[, cpos + seq.int(1L, m$S * ifelse(covariance, (m$S + 1L)/2L, 1L)) ],
                                     m$outputs, covariance)
   
   
   state = list(state.pred, state.var)
   names(state) <- c("sim", ifelse(covariance, "var", "sd"))
   
   output = list(output.pred, output.var)
   names(output) <- c("sim", ifelse(covariance, "var", "sd"))
   
   if (is.data.frame(newdata)) {
      pred <- list(state=state, output=output)      
   } else {
      pred <- vector(mode="list", length=length(newdata))
      Nprev <- 0L
      for (i in seq_along(newdata)) {
         idx <- seq.int(1L, nrow(newdata[[i]])) + Nprev
         Nprev <- Nprev + nrow(newdata[[i]])
         
         pred[[i]] <- list(state  = list(pred = state.pred[idx, , drop=FALSE],
                                         if (covariance) {
                                            state.var[, , idx]
                                         } else {
                                            state.var[idx, , drop=FALSE]
                                         }
         ),
         output = list(pred = output.pred[idx, , drop=FALSE],
                       if (covariance) {
                          output.var[, , idx]
                       } else {
                          output.var[idx, , drop=FALSE]
                       }
         )
         )
         names(pred[[i]]$state) <- c("sim",ifelse(covariance,"var","sd"))
         names(pred[[i]]$output) <- c("sim",ifelse(covariance,"var","sd"))
      }
      names(pred) <- names(newdata)
   }
   
   pred
}

ValidateResultObject <- function(object) {
   
   # Is there a CTSM model?
   if (!("model" %in% names(object) && class(object$model) == "ctsm"))
      stop("No CTSM model found in ", substitute(object))
   
   # Is xm there?
   if (!("xm" %in% names(object) && is.numeric(object$xm) && length(object$xm) == length(object$model$pars)))
      stop("No xm found in ", substitute(object))
}

PrepareEnv <- function(object, data, firstorderinputinterpolation, job, outputsize, covariance = 0, x0 = NULL, vx0 = NULL) {
   
   m <- object$model
   
   # Validate data
   data <- ValidateData(data, m$inputs, m$outputs, firstorderinputinterpolation)
   imat <- PackData(data, m$inputs, m$outputs)
   
   rho = new.env(parent=baseenv())
   
   m$options$sampletime <- as.numeric(!(attr(data, "ConstantSamplingTime")))
   m$options$inputinterpolation <- attr(data, "FirstOrderInputInterpolation")
   
   assign(".nobs", as.integer(attr(imat, "Observations")), rho)
   assign(".nset", as.integer(length(data)), rho)
   
   ints <- genInts(NULL, m, pred=TRUE, vx0 = !is.null(vx0))
   assign(".ints", as.integer(ints), rho)
   assign(".nint", as.integer(length(ints)), rho)
   
   PV.tmp <- data.frame(initial=object$xm)
   
   if (!is.null(x0)) {
      
      if (is.data.frame(x0)) { x0 <- unlist(x0) }
      if (is.matrix(x0)) { x0 <- x0[1,] }
      
      x0.names <- names(x0)
      if ( is.null(x0.names) || "" %in% x0.names ) {
         stop("x0 must be a named vector.", call. = FALSE)
      }
      
      x0.names[m$states %in% x0.names] <- paste0(x0.names[m$states %in% x0.names], "0")
      
      PV.tmp[x0.names, "initial"] <- x0
   }
   
   doubls <- genDoubles(PV.tmp, m, pred=TRUE)
   assign(".doubls", as.double(doubls), rho)
   assign(".ndoub", as.integer(length(doubls)), rho)
   
   assign(".tmat", as.double(attr(imat, "time")), rho)
   assign(".ntmat", Nt <- as.integer(sum(attr(imat, "Observations"))), rho)
   
   assign(".imat", as.double(imat), rho)
   assign(".nimat", as.integer(length(imat)), rho)
   
   omat <- numeric(Nt * outputsize)
   omat[1] <- covariance
   assign(".omat", omat, rho)
   assign(".nomat", as.integer(Nt * outputsize), rho)
   
   assign(".job", job, rho)
   
   assign(".info", 0L, rho)
   
   if (is.null(vx0)) {
      assign(".vx0", numeric(m$N*(m$N+1)/2), rho)
   } else {
      
      if (is.vector(vx0)) {
         if (length(vx0) != m$N*(m$N+1)/2) {
            stop("VX0 is not a covariance matrix")
         }
         
         vx0.matrix <- ToCovariance(m$N, matrix(vx0, 1))[, , 1]      
      }
      
      if (is.matrix(vx0)) {
         vx0.matrix <- vx0
         tryCatch(chol(vx0.matrix), 
                  error = function(e) {
                     stop("VX0 is not a covariance matrix", call.=FALSE)
                  }
         )
         
         vx0 <- vx0.matrix[lower.tri(vx0, diag = TRUE)]
      } else {
         stop("VX0 is not a covariance matrix")
      }
      
      assign(".vx0", vx0, rho)
   }
   
   return(rho)
}

CtsmrMethod <- function(object, job, outputsize, newdata, firstorderinputinterpolation) {
   ValidateResultObject(object)
   
   if (missing(newdata) || is.null(newdata)) {
      # Use data used while fitting
      newdata <- object$data
      constantsamplingtime <- attr(object$data, "ConstantSamplingTime")
      firstorderinputinterpolation <- attr(object$data, "FirstOrderInputInterpolation")
   }
   
   e <- PrepareEnv(object, newdata,
                   firstorderinputinterpolation = firstorderinputinterpolation,
                   job = job,
                   outputsize = outputsize)
   
   m <- object$model
   
   CallFunction("ctsmllike", m, e)
   
   if ((info <- get(".info", e)) != 0)
      stop(infocode(info), call.=FALSE, domain=NA)
   
   omat <- matrix(get(".omat", e), ncol=outputsize, byrow=TRUE)
   
   return( omat )
}