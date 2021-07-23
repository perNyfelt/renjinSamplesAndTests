library("codetools")

StandardAssignment <- setRefClass(
  Class ="StandardAssignment",
  fields = list(
    m_attributes = "list"
  ),
  methods = list(
    setAttribute = function(name, value) {
      m_attributes[[as.character(name)]] <<- as.character(value)
      return(invisible(.self))
    },
    getAttribute = function(name) {
      m_attributes[[name]]
    }
  )
)

SelfAssignment <- setRefClass(
  Class ="SelfAssignment",
  fields = list(
    m_attributes = "list"
  ),
  methods = list(
    setAttribute = function(name, value) {
      .self$m_attributes[[as.character(name)]] <- as.character(value)
      return(invisible(.self))
    },
    getAttribute = function(name) {
      m_attributes[[name]]
    }
  )
)

