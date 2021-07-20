library("se.alipsa:ctsmr")

model <- ctsm()
model$addSystem(dx ~ theta * (b - x) * dt + exp(sigma)*dw1)
model$addObs(y ~ x)
model$setVariance(yy ~ exp(S))
print(model)