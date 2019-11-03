CheckSystemEquation <- function(form, .self) {
   
   lhs <- form[[2]]
   rhs <- form[[3]]
   if (length(lhs) == 1 && is.name(lhs)) {
      # We have a single term on LHS
      state <- sub("^d?([[:alnum:]]+)", "\\1", lhs)
      states <- unique(c(.self$states, state))
   } else {
      stop("Incorrect equation. The left hand side must contain only a single term, eg. dXp")
   }
   
   #CheckExistance(state, .self)
   
   pars.all <- all.vars(rhs,unique=FALSE)
   pars <- unique(pars.all)
   diff.processes <- c("dt", pars[grep(pattern="d[(w\\d+)]",pars)])
   
   diff.terms <- lapply(diff.processes, function(x) { D(rhs, x) })
   names(diff.terms) <- diff.processes
   
   valid <- all(unlist(lapply(diff.terms,
                              function(x) { all(is.na(match(diff.processes, all.vars(x)))) })))
   
   if (!valid) { stop("Illegal cross terms between dt and dw") }
   
   # Find all parameters in the separated parts
   pars.after <- c(diff.processes, unlist(lapply(diff.terms, all.vars), use.names=FALSE))
   
   # Did we loose any terms like c in: f(.)*dt + g(.)*dw + c
   if (any(is.na(match(pars, pars.after)))) {
      stop("Illegal terms outside of the drift and diffusion terms")
   }
   
   # Illegal state dependence in the diffusion part g(.)
   valid <- all(unlist(lapply(diff.terms[!names(diff.terms) == "dt"], 
                              function(x) { all(is.na(match(states, all.vars(x)))) })))
   
   if (!valid) {
      stop("Illegal state dependence in the diffusion terms.")
   }
   
   
   res <- list(list(input = form, parsed = diff.terms, invariant = !any(pars == "t")))
   names(res) <- state
   return(res)
}

CheckObservationEquation <- function(form, .self) {
   lhs <- form[[2]]
   rhs <- form[[3]]
   
   if (!(length(lhs) == 1 && is.name(lhs))) {
      stop("Incorrect equation. The left hand side must contain only a single term, eg. y")
   }
   
   lhs <- as.character(lhs)
   
   # We have a single term on LHS
   
   #CheckExistance(lhs, .self)
   
   pars.all <- all.vars(rhs,unique=FALSE)
   pars <- unique(pars.all)
   
   illegal.pars <- pars[grep(pattern="d[(w\\d+)t]",pars)]
   if (!identical(pars[grep(pattern="d[(w\\d+)t]",pars)], character(0))) {
      stop("Illegal use of ", paste(illegal.pars, collapse=", "), " in the observation equation ", deparse(form), domain=NA)
   }
   
   res <- list(list(input=form, parsed=rhs, invariant=!any(pars == "t")))
   names(res) <- lhs
   return(res)
}

CheckExistance <- function(newvar, .self) {
   valid <- grepl(pattern="^[[:alpha:]][[:alnum:]]*$", x=newvar)
   if (!valid) { stop("Invalid variable name:", newvar) }
   
   if (any(newvar == .self$states)) { stop(newvar, " is already a state.", call.=FALSE, domain=NA) }
   
   if (any(newvar == .self$inputs)) { stop(newvar, " is already an input.", call.=FALSE, domain=NA) }
   
   if (any(newvar == .self$outputs)) { stop(newvar, " is already an output.", call.=FALSE, domain=NA ) }
}

ParseOutputVariance <- function(form, outputs) {
   
   lhs <- form[[2]]
   rhs <- form[[3]]
   
   ok <- FALSE
   
   output.vars <- lapply(outputs, as.name)
   N <- length(outputs)
   
   i <- 0
   while (!ok && i < N) {
      i <- i + 1
      j <- 0
      while (!ok && j < N) {
         j <- j + 1
         if (j == i) {
            combis <- c( output.vars[[i]],
                         call("^", output.vars[[i]], 2), 
                         call("*", output.vars[[i]],  output.vars[[i]]), 
                         as.name(paste(output.vars[[i]], output.vars[[i]], sep="")) )
         } else {
            combis <- c( call("*", output.vars[[i]], output.vars[[j]]), 
                         as.name(paste(output.vars[[i]], output.vars[[j]], sep="")) )
         }
         
         if (any(sapply(combis, identical, y=lhs))) {
            # We got a winner
            ok <- TRUE
         }
      }
   }
   
   if (ok) {
      return(list(i=i,j=j,var=rhs,input=form))
   }
   
   return(NULL)
}

ExpandEquations = function(eqs) {
   eqs$parsed <- eqs
   # Prepare the algebraic equations
   name_al <- names(eqs)
   alg_dep <- lapply(eqs, function(x) {name_al %in% all.vars(x)})
   notallowed <- NULL
   for (i in seq_along(name_al)) {
      # Did we already fix the equation?
      if (any(i == notallowed)) { next }
      # Check dependencies first
      eqs$parsed <- checkDependence(eqs$parsed, alg_dep, i)
   }
   return(eqs)
}

IsLinear <- function(form, .self) {
   states <- .self$states
   inputs <- .self$inputs
   
   A <- lapply(states, function(x) { stats::D(form, x)})
   names(A) <- states
   
   B <- lapply(inputs, function(x) { stats::D(form, x)})
   names(B) <- inputs
   
   # Check for linearity in states and inputs - don't care where the non-linearity is
   #    all(unlist(lapply(vars, function(x) { all(is.na(pmatch(vars, D(form, x))))})))
   vars <- c(states, inputs)
   # linear <- all(unlist(lapply(vars, function(x) { all(is.na(pmatch(vars, c(A, B))))})))
   
   linear <- all(is.na(match(vars, unlist(lapply(c(A,B),all.vars)))))
   
   # States and inputs present in form
   vars.present <- vars[!is.na(match(vars, all.vars(form)))]
   
   # pars.after <- unique(c(vars.present, unlist(lapply(c(A,B), all.vars), use.names=FALSE)))
   
   pars.after <- unique(c(vars.present, unlist(lapply(c(A,B), all.vars))))
   
   if (any(is.na(match(all.vars(form), pars.after)))) {
      # Parameters lost. Set to non-linear
      linear <- FALSE
   }
   
   if (!linear) { return(list(linear=linear)) }
   
   return(list(linear=linear, A=A, B=B))
}
