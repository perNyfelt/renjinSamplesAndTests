library('hamcrest')
library('se.alipsa:java-rc')

test.javaRc <- function() {
  greeter <- Greeter$new()
  greeter$setName("Per")
  assertThat(greeter$greet(), equalTo("Hello Per"))
}