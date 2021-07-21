# a pure R equivalent where the classes are implemented as S4 classes
# instead of java classes

employees <- list()
setClass("Employee", representation(name = "character", age = "numeric"))
employees <- append(employees, new("Employee", name = "Bob", age = 36))
employees <- append(employees, new("Employee", name = "Carol", age = 41))

age <- vector()
for (employee in employees) {
  age <- append(age, employee@age)
}
print(paste(
  "Average age is", mean(age), 
  "Number of employees are", length(age)
  )
)
