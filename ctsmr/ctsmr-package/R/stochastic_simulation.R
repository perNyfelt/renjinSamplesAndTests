#' Stochastic simulation for CTSM fits
#' @param fit Object of class 'ctsmr'
#' @param x0 x0
#' @param data data
#' @param dt Simulation dt for the brownian motion
#' @export
stochastic.simulate <- function(fit, x0, data, dt) {
   
   if (!is.data.frame(data)) {
      stop("data must be a data.frame")
   }
   
   m <- fit$model
   
   if (!all(m$inputs %in% colnames(data))) {
      stop("All inputs must be available in data")
   }
   
   time <- c("t","time")
   timei <- c("t","time") %in% colnames(data)
   
   if (all(timei)) {
      stop("The time must be called either t or time.")
   }
   
   t <- data[,time[timei]]
   
   uo <- c(t(as.matrix(data[,m$inputs])))
   
   nt <- nrow(data)
   
   x <- matrix(0, nrow=nt, ncol=m$N)
   y <- matrix(0, nrow=nt, ncol=m$S)
   
   out <- .Fortran("simulate", as.double(x), as.double(y), as.double(x0), as.double(uo), as.double(fit$xm),
                   as.double(t), as.integer(nt), as.double(dt),
                   PACKAGE = attr(m$libfile, "DLLname"))
   
   x <- as.data.frame(matrix(out[[1]], ncol=m$N, byrow=TRUE)[-nt,])
   colnames(x) <- m$states
   attr(x, "dt") <- dt
   
   y <- as.data.frame(matrix(out[[2]], ncol=m$S, byrow=TRUE)[-nt,])
   colnames(y) <- m$outputs
   
   list(state = x, output = y)
}