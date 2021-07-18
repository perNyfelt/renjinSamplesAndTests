library('hamcrest')
library('se.alipsa:studentS4')
library(stats)

test.creation <- function() {
  student1 <- Student("Per")
  result1 <- TestResult("Algebra", 78)
  assertThat(result1@grade, equalTo(78))

  # S4 objects are immutable so we need to reassign
  student1 <- addResult(student1, result1)
  assertThat(calculateGPA(student1), equalTo(78))
  
  student1 <- addResult(student1, TestResult("UX design", 88))
  assertThat(calculateGPA(student1), equalTo(83))
}
