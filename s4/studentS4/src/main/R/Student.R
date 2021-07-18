# remember to add export(function name) to NAMESPACE to make them available
Student <- setClass("Student",

  slots =  c(
    name = "character",
    results = "list"
  )

)

Student <- function(name, results = list()) {
  new("Student", name = name, results = results)
}

setGeneric("results", function(x) standardGeneric("results"))
setMethod("results", "Student", function(x) x@results)

setGeneric("addResult", function(x, value) standardGeneric("addResult"))
setMethod("addResult", "Student", function(x, value) {
  if (is.null(x@results)) {
    x@results <- c(value)
  } else {
    x@results <- c(x@results, value)
  }
  return(x)
})

setGeneric("calculateGPA", function(x) standardGeneric("calculateGPA"))
setMethod("calculateGPA", "Student", function(x) {
  tot <- 0
  count <- 0
  for (result in x@results) {
    tot <- tot + result@grade
    count <- count + 1 
  }
  return(tot/count)
})

setMethod("show", "Student", function(x) {
  cat("Student: ")
  cat(x@name)
  cat(", ")
  cat(as.character(length(x@results)))
  cat(" results\n")
})

setMethod("print", "Student", function(x) {
  cat("Student: ")
  cat(x@name)
  cat(", ")
  cat(as.character(length(x@results)))
  cat(" results\n")
})