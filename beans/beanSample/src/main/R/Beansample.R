# remember to add export(function name) to NAMESPACE to make them available
customerList <- list()

getCustomerList <- function() {
  return(customerList)
}

createCustomer <- function(name, age) {
  import(beans.Customer)
  Customer$new(name = name, age = age)  
}

addCustomer <- function(name, age) {
  # add to global var
  customer <- createCustomer(name, age)
  customerList <<- c(customerList, customer)
  return(customer)
}