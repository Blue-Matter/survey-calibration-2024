
# Simple example figure
m <- 1

x <- 1:20
xchange <- 10
x_fit <- x <= xchange
y <- m * x + ifelse(x <= xchange, 0, 7.5)

set.seed(42)
ysamp <- y + rnorm(length(y))

fit <- lm(ysamp ~ x + x_fit)
fit_notident <- lm(ysamp ~ factor(x))

png("manuscript/figures/simple_example.png", height = 2.5, width = 3.5, units = "in", res = 400)
par(mar = c(2, 2, 1, 1), mfrow = c(1, 1))
plot(x, ysamp, xlab = NULL, ylab = NULL, pch = 16, panel.first = grid(),
     xaxt = "n", yaxt = "n")
axis(1, labels = FALSE)
axis(2, labels = FALSE)
mtext("Year", 1, line = 1)
mtext("Index", 2, line = 1)

#lines(x, predict(fit_notident))

lines(x[x_fit], predict(fit)[x_fit])
lines(x[!x_fit], predict(fit)[!x_fit])

#lines(x[x_fit], y[x_fit])
#lines(x[!x_fit], y[!x_fit])
abline(v = xchange + 0.5, lty = 2)
dev.off()


png("manuscript/figures/simple_example2.png", height = 2.5, width = 6, units = "in", res = 400)
par(mar = c(2, 2, 2, 1), mfrow = c(1, 2))

plot(x, ysamp, xlab = NULL, ylab = NULL, pch = 16, panel.first = grid(),
     xaxt = "n", yaxt = "n")
axis(1, labels = FALSE)
axis(2, labels = FALSE)
mtext("Year", 1, line = 1)
mtext("Index", 2, line = 1)
title(expression(rho~~"is identifiable"))

lines(x[1:xchange], y[1:xchange])
lines(x[-c(1:xchange)], y[-c(1:xchange)])
abline(v = xchange + 0.5, lty = 2)

plot(x, ysamp, xlab = NULL, ylab = NULL, pch = 16, panel.first = grid(),
     xaxt = "n", yaxt = "n")
axis(1, labels = FALSE)
axis(2, labels = FALSE)
mtext("Year", 1, line = 1)
mtext("Index", 2, line = 1)

title(expression(rho~~"is not identifiable"))

lines(x, predict(fit_notident))

abline(v = xchange + 0.5, lty = 2)

dev.off()

