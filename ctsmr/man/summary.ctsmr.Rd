\name{summary.ctsmr}
\alias{summary.ctsmr}
\title{Simple summary of a CTSMR class}
\usage{
\method{summary}{ctsmr}(object, correlation = extended, extended = FALSE,
  ...)
}
\arguments{
  \item{object}{an object of class "ctsmr".}

  \item{correlation}{logical. if TRUE, the correlation
  matrix of the estimated parameters is returned and
  printed.}

  \item{extended}{logical. If TRUE, additional information
  is returned and printed. Implies correlation = TRUE.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\description{
Simple summary of a CTSMR class
}

