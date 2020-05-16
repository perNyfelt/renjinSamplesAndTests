Greeter <- setRefClass("Greeter",
  
  fields = list (
    hello = "ANY"
  ),
   
  methods = list(
    
    initialize = function() {     
      hello <<- Hello$new()
    },
    
    setName = function(name) {
      hello$setName(name)
    },
    
    greet = function() {
      hello$sayHello()
    } 
  )
)