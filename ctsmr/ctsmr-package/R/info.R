
infocode <- function(info) {
   
   info <- as.integer(info)
   
   if (info == 0)
      msg <- "Converged."
   else if (info == -1)
      msg <- "Terminated"
   else if (info == 2)
      msg <- "The maximum number of objective function evaluations has been exceeded."
   else if (info == 5)
      msg <- "The prior covariance matrix is not positive definite."
   else if (info == 10)
      msg <- "The amount of data available is insufficient to perform the estimation."
   else if (info == 20)
      msg <- "The maximum objective function value (1e300) has been exceeded."
   else if (info == 30)
      msg <- "The state covariance matrix is not positive definite."
   else if (info == 40)
      msg <- "The measurement noise covariance matrix is not positive definite."
   else if (info == 50)
      msg <- "Unable to calculate matrix exponential."
   else if (info == 60)
      msg <- "Unable to determine reciprocal condition number."
   else if (info == 70)
      msg <- "Unable to compute singular value decomposition."
   else if (info == 80)
      msg <- "Unable to solve system of linear equations."
   else if (info == 90)
      msg <- "Unable to perform numerical ODE solution."
   else
      msg <- "I'm sorry, but I don't know this error code."

   
   msg <- paste(msg, "Code:", info)
   return(msg)
}