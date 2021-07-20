
setOldClass("DLLInfo")


# Options and default values

.defaultoptions <- list(
   # Filter options
         solutionMethod         = 1,     # Method for solution propagation: subsampling, Adams, BDF
         InitialVarianceScaledIdentity = FALSE,
         initialVarianceScaling = 1.0,     #
         numberOfSubsamples     = 10,      # Valid for LTI and NL
         odeeps                 = 1.0e-12, #
         nIEKF                  = 1,       #
         iEKFeps                = 1.0e-12, #
   # Optimization options
         maxNumberOfEval        = 500,     # Maximum number of objective function evaluations
         eps                    = 1.0e-6,  # Adjustment factor for initial step length in line search
         eta                    = 1.0e-14,  # Relative error in calculation of objective function
   # Computational options
         hubersPsiLimit         = 3.0,     # Cut-off value for Huber's psi-function
         padeApproximationOrder = 6,       # Pade approximation order
         svdEps                 = 1.0e-12, # Tolerance for singular value decomposition
   # Advanced optimization options
         lambda                 = 1.0e-4,  # Lagrange multiplier in penalty function
         smallestAbsValueForNormalizing = 1.0e-30 # Minimum absolute value used for normalizing in penalty function
      )

CTSM_NA_VALUE <- 1e+300