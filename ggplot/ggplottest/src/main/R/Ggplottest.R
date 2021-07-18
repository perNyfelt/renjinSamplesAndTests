library(ggplot2)
png()

ggplot(ToothGrowth, aes(supp, len)) + geom_violin()
dev.off()