library('hamcrest')

test.selfAssignment <- function() {
  s1 <- StandardAssignment$new()
  s1$setAttribute("name", "foo")
  assertThat(s1$getAttribute("name"), equalTo("foo"))
}

test.standardAssignment <- function() {
  s2 <- SelfAssignment$new()
  s2$setAttribute("name", "foo")
  assertThat(s2$getAttribute("name"), equalTo("foo"))
}


