library('se.alipsa:beanSample')

# as long as we are in the same session the customer list will keep increasing
# Reset it to give us a fresh start
resetCustomerList()
addCustomer(name = "Bob", age = 36)
addCustomer(name = "Carol", age = 41)

age <- vector()
for (customer in getCustomerList()) {
  age <- append(age, customer$age)
}
print(paste(
  "Average age is", mean(age), 
  "Number of customers are", length(age)
  )
)

