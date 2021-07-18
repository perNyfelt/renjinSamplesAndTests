
TestResult <- setClass("TestResult",

  slots = c(
    subject = "character",
    grade = "numeric"
  ),
  prototype = list(
    subject = NA_character_,
    grade = NA_real_
  )
)

TestResult <- function(subject = NA, grade = 0) {
  new("TestResult", subject = subject, grade = grade)
}

setMethod("show", "TestResult", function(x) {
  cat("TestResult:")
  cat(x@subject)
  cat(", ")
  cat(as.character(x@grade))
  cat("\n")
})

setMethod("print", "TestResult", function(x) {
  cat("TestResult:")
  cat(x@subject)
  cat(", ")
  cat(as.character(x@grade))
  cat("\n")
})