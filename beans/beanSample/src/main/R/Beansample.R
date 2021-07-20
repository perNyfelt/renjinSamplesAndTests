#
customerList <- list()

getCustomerList <- function() {
  return(customerList)
}

createCustomer <- function(name, age) {
  import(beans.Customer)
  Customer$new(name = name, age = age)  
}

addCustomer <- function(name, age) {
  customer <- createCustomer(name, age)
  # add to global list
  customerList <<- c(customerList, customer)
  return(customer)
}