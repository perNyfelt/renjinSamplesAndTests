
#' Simple summary of a CTSMR class
#' @method summary ctsmr
#' @S3method summary ctsmr
#' @param object an object of class "ctsmr".
#' @param correlation logical. if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param extended logical. If TRUE, additional information is returned and printed. Implies correlation = TRUE.
#' @param ... further arguments passed to or from other methods.
summary.ctsmr <- function(object, correlation=extended, extended = FALSE, ...) {
   
   stopifnot(is.logical(correlation))
   stopifnot(is.logical(extended))
   
   z <- object
   
   ans <- list()
   #ans$coefficients <- cbind(z$xm,z$sd,z$t,z$pt)
   c <- matrix(NA,nrow=length(z$xm),ncol=ifelse(!extended, 4, 6))
   c[,1] <- z$xm
   c[z$estimated,2:4] <- cbind(z$sd,z$t,z$pt)
   
   if (extended)
      c[z$estimated, c(5,6)] <- cbind(z$dL, z$dpen)

   rownames(c) <- names(z$xm)
   colnames(c) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", 
                    if (extended)
                       c("dF/dPar", "dPen/dPar")
                    else
                       NULL
                    )
   
   ans$coefficients <- c
   
   if (is.logical(correlation) && correlation) {
      ans$correlation <- z$corr
   }
   
   class(ans) <- "summary.ctsmr"
   
   ans
}

#' Print the summary of a CTSMR class
#' @method print summary.ctsmr
#' @S3method print summary.ctsmr
#' @param x an object of class "summary.ctsmr".
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
print.summary.ctsmr <- function(x, digits = max(3, getOption("digits") - 3),...) {
   
   cat("Coefficients:\n")
   
   printCoefmat(x$coefficients)
   
   if (!is.null(x$corr)) {
      cat("\nCorrelation of coefficients:\n")
      
      # Do the same thing as stats:::print.summary.lm
      correl <- format(round(x$corr, 2), nsmall = 2, digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1, -(dim(x$corr)[1]), drop = FALSE], quote = FALSE)
   }
   
   cat("\n")
   invisible(x)
}

#' Print the CTSMR class
#' @method print ctsmr
#' @S3method print ctsmr
#' @param x an object of class "ctsmr".
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
print.ctsmr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
   cat("\nCTSM model:\n\n")
   x$model$show()
   
   cat("\nCoefficients:\n")
   
   print.default(format(x$xm, digits = digits), print.gap = 2, 
                 quote = FALSE)
   
   cat("\n")
   invisible(x)
}