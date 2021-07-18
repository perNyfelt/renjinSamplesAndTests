import(org.jfxplot.PlotApp)
import(org.jfxplot.GraphicsState)
import(org.jfxplot.StageManager)
import(org.jfxplot.PlotManager)
library('se.alipsa.jfxplot:plot')

x = seq(-2.0*pi,2.0*pi,4.0*pi/100)
y = sin(x)

data <- c(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0)

m <- matrix(data,3, 3)
m <- m * 100.0
rownames(m) <- c("A","B","C")
colnames(m) <- c("X","Y","Z")

plotManager <- PlotManager$new
stageManager <- StageManager$new
plot(x,y)
#dev.new()
#barplot(m)
#dev.new()
#canvas()