import(se.alipsa.charsequence.Concatenator)

concatenator <- Concatenator$new()
msg <- concatenator$concatString("hello ", "world")
print(msg)
msg <- concatenator$concatCharSeq("hello ", "world2")
print(msg)