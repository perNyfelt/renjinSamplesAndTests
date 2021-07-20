\name{filter.ctsmr}
\alias{filter.ctsmr}
\title{Optimal filter method for CTSM fits}
\usage{
filter.ctsmr(object, newdata = NULL, firstorderinputinterpolation = FALSE,
  covariance = FALSE, ...)
}
\arguments{
  \item{object}{Object of class 'ctsmr'}

  \item{newdata}{An optional data frame with new data.}

  \item{firstorderinputinterpolation}{If TRUE (default
  FALSE) first order linear interpolation of the inputs
  between samples.}

  \item{covariance}{If TRUE (default FALSE) the full
  covariance matrix will be returned.}

  \item{...}{Unused.}
}
\description{
Optimal filter method for CTSM fits
}

