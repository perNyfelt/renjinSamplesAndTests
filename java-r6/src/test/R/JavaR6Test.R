library('hamcrest')
library('se.alipsa:java-r6')

test.javaR6 <- function() {
  greeter <- Greeter$new()
  greeter$setName("Per")
  assertThat(greeter$greet(), equalTo("Hello Per"))
}