library(R6)

Greeter <- R6Class("Greeter",
  
  private = list (
    hello = "ExternalPtr"
  ),
   
  public = list(
    
    initialize = function() {     
      private$hello <- Hello$new()
    },
    
    setName = function(name) {
      private$hello$setName(name)
    },
    
    greet = function() {
      private$hello$sayHello()
    } 
  )
)