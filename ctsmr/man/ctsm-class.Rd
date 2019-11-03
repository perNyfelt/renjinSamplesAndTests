\docType{class}
\name{ctsm-class}
\alias{ctsm}
\alias{ctsm-class}
\title{CTSM class}
\description{
CTSM Class (description)
}
\details{
\code{ctsm} class is a reference class (details).
}
\section{Methods}{
\itemize{
\item \code{AnalyseModel()}: Analyses the structure of the model.

\item \code{addInput(...)}: Specifies which variables are inputs.

\item \code{addObs(form)}: Adds one or more observation equations to the model.

\item \code{addSystem(form)}: Adds one or more stochastic differential equations to the model. The diffusion processes must be called dw{1,2,...}

\item \code{conditionallikelihood(fit, which, maxpts = 100, data)}: Calculates a conditional likelihood.

\item \code{estimate(data, firstorderinputinterpolation = FALSE, threads = 1, debug = FALSE, verbose = FALSE)}: Estimate the model. Use threads >= 1 for parallel computation.

\item \code{pred(res, newdata, opt, k = NULL, covariance = FALSE)}: General prediction, smoothing, filtering method.

\item \code{profilelikelihood(fit, which, data, maxpts = 100, nsd = 5, threads = 1)}: Calculates a profile likelihood.

\item \code{setOptions(con = list(), default = FALSE)}: Use to change the default options.

\item \code{setParameter(...)}: Use to set the initial value and lower and upper bounds.

\item \code{setVariance(form)}: Sets the variance or covariance of the outputs.
}}
\section{Fields}{
  \itemize{ \item N. Number of states \item M. Number of
  inputs \item S. Number of outputs \item NPARAM. Number of
  parameters }
}

\section{Methods}{
  \itemize{ \item setOptions \item compile }
}

