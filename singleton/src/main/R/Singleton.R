import(se.alipsa.singleton.Validator)

validator <- Validator$getInstance()

validateString <- function(str) {
  validator$validate(str)
}  