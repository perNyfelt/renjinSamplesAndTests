library('hamcrest')
library('se.alipsa:beanSample')

test.beanSample <- function() {
  addCustomer(name = "Bob", age = 36)
  addCustomer(name = "Carol", age = 41)
  assertThat(length(getCustomerList()), equalTo(2))
  
  assertThat(getCustomerList()[[1]]$name, equalTo("Bob"))
  assertThat(getCustomerList()[[2]]$age, equalTo(41))
}

test.getterStillWorks <- function() {
  # We can create and add stuff to a java Map as an alternative to using R list()
  import(java.util.HashMap)
  ageMap <- HashMap$new()
  bobby <- createCustomer(name = "Bobby", age = 26)
  ageMap$put("bobby", bobby)
  bobby <- ageMap$get("bobby")
  assertThat(bobby$getName(), equalTo("Bobby"))
  assertThat(bobby$getAge(), equalTo(26))
}

test.explicitCreation <- function() {
  import(beans.Customer)
  bob <- Customer$new(name = "Bob", age = 36)
  carol <- Customer$new(name = "Carol", age = 41)
  cat("'bob' is an ", typeof(bob), ", bob$name = ", bob$name, "\n", sep = "")
  cat("'bob' is a ", class(bob), ", bob$getName() = ", bob$getName(), "\n", sep = "")
}