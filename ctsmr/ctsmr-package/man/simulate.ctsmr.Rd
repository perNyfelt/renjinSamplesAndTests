\name{simulate.ctsmr}
\alias{simulate.ctsmr}
\title{Simulate method for CTSM fits}
\usage{
\method{simulate}{ctsmr}(object, nsim = 1L, seed = NULL, newdata = NULL,
  firstorderinputinterpolation = FALSE, ...)
}
\arguments{
  \item{object}{Object of class 'ctsmr'}

  \item{nsim}{number of realisations. Defaults to 1.}

  \item{seed}{seed}

  \item{newdata}{An optional data frame with new data.}

  \item{firstorderinputinterpolation}{If TRUE (default
  FALSE) first order linear interpolation of the inputs
  between samples.}

  \item{...}{Unused.}
}
\description{
Simulate method for CTSM fits
}

