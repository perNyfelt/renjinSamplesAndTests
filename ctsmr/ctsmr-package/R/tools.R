genDoubles <- function(ParameterValues, m, pred=FALSE) {
   
   N <- m$N
   nparam <- m$NPARAM
   options <- m$options
   linear <- m$linear
   pars <- m$pars
   
   doubls <- numeric((2*nparam+9)*nparam+14)
   
   doubls[(nparam+5)*nparam+9]  <- options$svdEps
   doubls[(nparam+5)*nparam+11] <- options$initialVarianceScaling
   
   
   if (!pred) {
      doubls[1]                    <- options$eps
      doubls[2]                    <- options$eta
      doubls[(nparam+5)*nparam+7]  <- options$hubersPsiLimit
      doubls[(nparam+5)*nparam+8]  <- options$lambda
      
      doubls[(nparam+5)*nparam+10] <- options$smallestAbsValueForNormalizing
   }

   ## Prior correlation matrix
   ## doubls[(nparam+5)*nparam+11 +  (1:(npr*(npr+1)))] <- 
   
#    ### Initial state - i = 1:N
#    #i <- 1:.self$N
#    i <- match(estState$name,.self$pars[1:.self$N])
#    # Initial value
#    doubls[(2*nparam+6)*nparam+12+i] <- estState$initval
#    
#    if (!pred) {
#       # Upper bound
#       doubls[(2*nparam+7)*nparam+12+i] <- estState$ub
#       # Lower bound
#       doubls[(2*nparam+8)*nparam+12+i] <- estState$lb
#       # MAP - prior standard deviance
#       # j goes from 1 to #MAP vars: N_MAP_INITIAL
#       # j <- 1:N_MAP_INITIAL
#       # doubls[(2*nparam+5)*nparam+12+j] <-
#       # 0=fixed, 1=ML, 2=MAP
#    }
#    
#    ## Parameters - i = 1:11
#    #i <- 1:(nparam-n)
#    i <- match(estPar$name,.self$pars[-i])
#    # Initial value
#    doubls[(2*nparam+6)*nparam+12+n+i] <- estPar$initval
#    
#    if (!pred) {
#       # Upper bound
#       doubls[(2*nparam+7)*nparam+12+n+i] <- estPar$ub
#       # Lower bound
#       doubls[(2*nparam+8)*nparam+12+n+i] <- estPar$lb
#       # MAP - prior standard deviance
#       # j continues from N_MAP_INITIAL to N_MAP_PARAMETER
#       # j <- N_MAP_INITIAL+(1:N_MAP_PARAMETER)
#       # doubls[(2*nparam+5)*nparam+12+j] <-
#    }
   
   i <- 1:nparam
   doubls[(2*nparam+6)*nparam+12+i] <- ParameterValues[pars, "initial"]
   
   if (!pred) {
      # Upper bound
      doubls[(2*nparam+7)*nparam+12+i] <- ParameterValues[pars, "upper"]
      # Lower bound
      doubls[(2*nparam+8)*nparam+12+i] <- ParameterValues[pars, "lower"]
      doubls[is.na(doubls)] <- 0
      
      if (m$nprior > 0) {
         ## MAP
         estType <- ParameterValues[pars, "type"]
         priornames <- (pars)[match(estType, "map", 0L) == 1L]
         
         corr <- CheckCorrelationMatrix(m$PriorCorrelationMatrix, priornames)
         
         # Prior correlation matrix
         doubls[(nparam+5L)*nparam+11L+seq(1L,m$nprior^2)] <- c(corr)
         doubls[(2L*nparam+5L)*nparam+12L+seq(1L,m$nprior)] <- ParameterValues[priornames, "prior"]
      }
   }
   
   doubls[length(doubls)+(-1:0)] <- c(options$iEKFeps, options$odeeps)
   
   return (as.double(doubls))
}

genInts <- function(ParameterValues, m, pred=FALSE, vx0 = FALSE) {
   
   N <- m$N
   nparam <- m$NPARAM
   options <- m$options
   linear <- m$linear
   pars <- m$pars
   
   ints  <- integer(nparam+11L)
   
   if (!pred) {
      ints[1L] <- options$maxNumberOfEval
      
      estType <- ParameterValues[pars, "type"]
      est <- numeric(nparam)
      
      # Number of priors
      m$nprior <- sum(match(estType, "map", 0L))
      ints[9L] <- m$nprior

      # 0=fixed, 1=ML, 2=MAP
      est[estType == "fixed"] <- 0L
      est[estType == "ml"]    <- 1L
      est[estType == "map"]   <- 2L
      
      m$Nest <- sum(est != 0L)
      
      ints[10L + seq(1L,nparam)] <- est
   }
   
   # Constant sampling = 0. Variable sampling = 1
   ints[4L] <- options$sampletime
   # Hold order. Zero order = 0. First order = 1
   ints[5L] <- options$inputinterpolation
   

   # Solution method:
   #  - for NL
   # Subsampling approximation                                   : solutionMethod <- 0
   # Numerical ODE solution (Adams method for non-stiff systems) : solutionMethod <- 1
   # Numerical ODE solution (BDF for stiff systems)              : solutionMethod <- 2
   # - otherwise solutionMethod <- 0
   #
   # Type:
   #
   # LTI : type <- 0
   # LTV : type <- 1
   # NL  : type <- 2
   #
   ints[6] <- ifelse(linear == TRUE, 0, options$type + options$solutionMethod)
   # Padee approximation order
   ints[7] <- options$padeApproximationOrder
   # For subsampling approximation
   ints[8] <- ifelse(options$solutionMethod==0, options$numberOfSubsamples, 1)
   # If TRUE the initial variance matrix will be k*I
   ints[10L] <- ifelse(vx0, 2L, as.integer(options$InitialVarianceScaledIdentity))

   ints[length(ints)]            <- options$nIEKF

   return (ints)
}

codeWalker <- function(ex,node,...) {
   if (is.list(ex) && length(ex) > 0) {
      for (j in 1:length(ex)) 
         ex[[j]] <- Recall(ex[[j]],node,...)
   }
   if (typeof(ex) == "language") {
      
      exx <- as.list(ex)
      
      # Check if  [...] / -x
      if (length(exx) == 3L) {
         if (any(exx[[1L]] == c('/','*','^','+','-'))) {
            if (is.numeric(exx[[3L]]) && exx[[3L]] < 0) {
               exx[[3L]] <- as.call(list(as.symbol("("), exx[[3L]]))
            }
         }
      }
      
      for (i in 1:length(exx)) {
         ex[[i]] <- Recall(exx[[i]],node,...)
      }
   }
   # Reached a node
   ex <- node(ex,...)
   return (ex)
}

insertEquations <- function(ex,eqs) {
   if (typeof(ex) == "symbol") {
      exc <- as.character(ex)
      if (exc %in% names(eqs))
         ex <- eqs[[exc]]
   }
   return (ex)
}

convertToFortranSyntax <- function(ex,xm,xp,uo) {
   if (typeof(ex) == "symbol") {
      exc <- as.character(ex)
      if (exc %in% xm)
         return (as.call(list(as.name("XM"),as.numeric(which(exc==xm)))))
      if (exc %in% xp)
         return (as.call(list(as.name("XP"),as.numeric(which(exc==xp)))))
      if (exc %in% uo)
         return (as.call(list(as.name("UO"),as.numeric(which(exc==uo)))))
      if (exc %in% c("t", "T"))
         return (as.symbol("T"))
      #if (exc %in% names(eqs))
      #   return (eqs[[exc]])
   }
   if (is.numeric(ex)) {
      #return (as.name(gsub("\\.?0*E","D",sprintf("%.14E",ex))))
      return (gsub("e","D",format(ex,digits=14,scientific=TRUE)))
   } else if (is.null(ex)) {
      return ("0D+00")
   }
   
   return (ex)
}


diffVectorFun <- function(f,x) {
   nf <- length(f)
   nx <- length(x)
   # Output as a list
   jac <- vector("list",nf*nx)
   
   # Column major
   k <- 0
   for (j in 1:nx)
      for (i in 1:nf)
         jac[[k<-k+1]] <- stats::D(f[[i]],x[j])
   jac
}

getPars <- function(eq) {
   if (is.list(eq))
      unique(unlist(lapply(eq,all.vars)))
   else
      unique(all.vars(eq))
}

convertToFortran <- function(lfcn,aux,N,M,append,pars,states,uovars) {
   code <- aux$header
   k <- length(code)
   lff <- codeWalker(lfcn,convertToFortranSyntax,pars,states,uovars)
   # Remove spaces and "s introduced by convertToFortranSyntax
   lff <- lapply(lff, function(x){gsub("\\^","**",gsub('[[:space:]"]',"",paste(deparse(x),collapse="")))})
   for (j in seq_len(M)) {
      for (i in seq_len(N)) {
         code[[k<-k+1]] <- sprintf("%s(%i,%i)=%s",aux$matname,i,j,
                       if (append==1 && j == M) { "0D+00" } else { lff[[(j-1)*N+i]] })
      }
   }
   if (append == 0) {
      for (i in seq_len(N))
         code[[k<-k+1]] <- sprintf("%s(%i)=%s",aux$funname,i,lff[[M*N+i]])
   }
   code[[k+1]] <- "END"
   return (code)
}

genGlobal <- function(N,M,S,NPARAM) {
   code <- list("INTEGER N,M,S,NPARAM",
                sprintf("PARAMETER (N=%i)",N),
                sprintf("PARAMETER (M=%i)",M),
                sprintf("PARAMETER (S=%i)",S),
                sprintf("PARAMETER (NPARAM=%i)",NPARAM))
   return (code)
}

addNames <- function(x, name) {

   if (is.matrix(x)) {
      if (dim(x)[1] != length(name)) {
         warning(paste("Internal error! Number of names and numbers doesn't match. Dimension of the matrix = ", dim(x)[1], "x", dim(x)[2]," and number of names = ", length(name), sep=""))
         return(x)
      }
      colnames(x) <- rownames(x) <- name
   } else {
      if (length(x) != length(name)) {
         warning(paste("Internal error! Number of names and numbers doesn't match. Length = ", length(x), " and number of names = ", length(name), sep=""))
         return(x)
      }
      names(x) <- name
   }
   
   return(x)
}


.fvecxj.aux<-list(header=list("SUBROUTINE FVECXJ(XM,XP,UO,T,JACOB,F,NPARAM,N,M)",
                              "INTEGER NPARAM,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,JACOB(N,N),F(N)"),
                  matname="JACOB",
                  funname="F")

.fvecuj.aux<-list(header=list("SUBROUTINE FVECUJ(XM,XP,UO,T,JACOB,NPARAM,N,M)",
                              "INTEGER NPARAM,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,JACOB(N,M)"),
                  matname="JACOB")

.hvecxj.aux<-list(header=list("SUBROUTINE HVECXJ(XM,XP,UO,T,JACOB,H,NPARAM,N,M,S)",
                              "INTEGER NPARAM,N,M,S",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,JACOB(S,N),H(S)"),
                  matname="JACOB",
                  funname="H")

.hvecuj.aux<-list(header=list("SUBROUTINE HVECUJ(XM,XP,UO,T,JACOB,NPARAM,N,M,S)",
                              "INTEGER NPARAM,N,M,S",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,JACOB(S,M)"),
                  matname="JACOB")

.smat.aux  <-list(header=list("SUBROUTINE SMAT(XM,UO,T,SM,NPARAM,N,M,S)",
                              "INTEGER NPARAM,S,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,SM(S,S)"),
                  matname="SM")

.sigmat.aux<-list(header=list("SUBROUTINE SIGMAT(XM,UO,T,SIGMA,NPARAM,N,M)",
                              "INTEGER NPARAM,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,SIGMA(N,N)"),
                  matname="SIGMA")

.amat.aux <- list(header=list("SUBROUTINE AMAT(XM,XP,UO,T,A,NPARAM,N,M)",
                              "INTEGER NPARAM,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,A(N,N)"),
                  matname="A")

.bmat.aux <- list(header=list("SUBROUTINE BMAT(XM,XP,UO,T,B,NPARAM,N,M)",
                              "INTEGER NPARAM,N,M",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,B(N,M)"),
                  matname="B")

.cmat.aux <- list(header=list("SUBROUTINE CMAT(XM,XP,UO,T,C,NPARAM,N,M,S)",
                              "INTEGER NPARAM,N,M,S",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,C(S,N)"),
                  matname="C")

.dmat.aux <- list(header=list("SUBROUTINE DMAT(XM,XP,UO,T,D,NPARAM,N,M,S)",
                              "INTEGER NPARAM,N,M,S",
                              "DOUBLE PRECISION XM(NPARAM),XP(N),UO(M),T,D(S,M)"),
                  matname="D")
