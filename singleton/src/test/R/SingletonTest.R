library('hamcrest')


print(paste("123 is valid =", validateString("123")))

assertThat(validateString("123"), equalTo(TRUE))