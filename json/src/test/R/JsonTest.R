library('hamcrest')

saveRDS(c(1,2,3), "robj.rds")
  

test.rjson <- function() {
  robj <- readRDS("robj.rds")
  
  library('org.renjin.cran:rjson')
  jobj <- toJSON(robj)
  print(jobj)
  #unloadNamespace("rjson") # not implemented in renjin
  detach("package:rjson", unload=TRUE) # only half working in renjin
}

test.jsonlite <- function() {
  robj <- readRDS("robj.rds")
  library('org.renjin.cran:jsonlite')
  jobj <- toJSON(robj)
  print(jobj)
  detach("package:jsonlite", unload=TRUE) # only half working in renjin
}

# NOT WORKING
#test.jsonify <- function() {
#  robj <- readRDS("robj.rds")
#  library('org.renjin.cran:jsonify')
#  jobj <- to_json(robj)
#  print(jobj)
#  detach("package:jsonify", unload=TRUE) # only half working in renjin
#}