\name{predict.ctsmr}
\alias{predict.ctsmr}
\title{Predict method for CTSM fits}
\usage{
\method{predict}{ctsmr}(object, n.ahead = 1, covariance = FALSE,
  newdata = NULL, firstorderinputinterpolation = FALSE, x0 = NULL,
  vx0 = NULL, ...)
}
\arguments{
  \item{object}{Object of class 'ctsmr'}

  \item{n.ahead}{The number of steps ahead for which
  prediction is required.}

  \item{covariance}{Should the full covariance matrix be
  outputted.}

  \item{newdata}{An optional data frame with new data.}

  \item{firstorderinputinterpolation}{If TRUE (default
  FALSE) first order linear interpolation of the inputs
  between samples.}

  \item{x0}{initial state}

  \item{vx0}{initial state covariance}

  \item{...}{Unused.}
}
\description{
Predict method for CTSM fits
}

