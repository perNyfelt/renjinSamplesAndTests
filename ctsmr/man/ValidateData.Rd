\name{ValidateData}
\alias{ValidateData}
\title{Validates user given data for CTSM}
\usage{
ValidateData(data, inputs, outputs, firstorderinputinterpolation)
}
\arguments{
  \item{data}{data.frame or list of data.frames.}

  \item{inputs}{vector of names of the input variables.}

  \item{outputs}{vector of names of the output variables.}

  \item{firstorderinputinterpolation}{logical. If TRUE, a
  first order interpolation between input samples will be
  used. Otherwise constant hold.}
}
\description{
Validates user given data for CTSM
}

