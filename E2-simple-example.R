
# Simple example figure
m <- 1

x <- 1:20
xchange <- 10
y <- m * x + ifelse(x <= xchange, 0, 7.5)

set.seed(42)
ysamp <- y + rnorm(length(y))

png("manuscript/figures/simple_example.png", height = 2.5, width = 3.5, units = "in", res = 400)
par(mar = c(2, 2, 1, 1))
plot(x, ysamp, xlab = NULL, ylab = NULL, pch = 16, panel.first = grid(),
     xaxt = "n", yaxt = "n")
axis(1, labels = FALSE)
axis(2, labels = FALSE)
mtext("Year", 1, line = 1)
mtext("Index", 2, line = 1)


lines(x[1:xchange], y[1:xchange])
lines(x[-c(1:xchange)], y[-c(1:xchange)])
abline(v = xchange + 0.5, lty = 2)
dev.off()
