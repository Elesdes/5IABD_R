#1) Analyse de Fourier
# -> Série de fourrier (SF)
# -> Transformation de fourier continue (TF)
# -> Transformation de fourier discrète (TFD)+ TFD Rapide


f <- function(t) {
 
  return((pi  - abs(t %% (2 * pi))))   
}

f <- function(t) { ifelse(t >= -pi & t <= pi, pi - abs(t), pi - abs(t - 2*pi*floor((t + pi)/(2*pi)))) }


g <- function(t) { ifelse(t >= -pi & t <= pi, t^2 - pi^2, (t - 2*pi*floor((t + pi)/(2*pi)))^2 - pi^2) }

h <- function(t) {
  if (t >= 0 && t < pi) {
    return(exp(-t / pi))
  } else if (t >= pi && t < 2 * pi) {
    return(0)
  }
}

h <- function(t) {
  t <- t %% (2*pi)
  ifelse(t >= 0 & t < pi, exp((-t )/pi), 0) 
}
x <- seq(-5 * pi, 5 * pi, by= 0.01)

plot(x, f(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")


plot(x, g(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")

plot(x, h(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "green")
#-----------------------------------------------------------------------------------------------------------------------------------
#2)

# Calculate the Fourier series coefficients for f

N <- as.numeric(readline("Saisissez la valeur de N :"))
N <- 5
a0 <- 1/(2*pi) * integrate(f, -pi, pi)$value
an <- rep(0, N)
bn <- rep(0, N)
iOne <- rep(0, N)
iTwo <- rep(0, N)
cn <- rep(0, N)

for (n in 1:N) {
  an[n] <- 1/pi * integrate(function(t) f(t) * cos(n*t), -pi, pi)$value
  bn[n] <- 1/pi * integrate(function(t) f(t) * sin(n*t), -pi, pi)$value
  iOne[n] <- integrate(function(t) f(t) * cos(n*t), -pi, pi)$value
  iTwo[n] <- integrate(function(t) f(t) * sin(n*t), -pi, pi)$value
  cn[n] <- iOne[n] + i * iTwo[n]
}
a0
an[N]
bn[N]
cn[N]
# Calculate the Fourier series coefficients for g
for (n in 1:N) {
  an[n] <- 1/pi * integrate(function(t) g(t) * cos(n*t), -pi, pi)$value
  bn[n] <- 1/pi * integrate(function(t) g(t) * sin(n*t), -pi, pi)$value
  iOne[n] <- integrate(function(t) g(t) * cos(n*t), -pi, pi)$value
  iTwo[n] <- integrate(function(t) g(t) * sin(n*t), -pi, pi)$value
  cn[n] <- iOne[n] + i * iTwo[n]
}
a0
an[N]
bn[N]
cn[N]

# Calculate the Fourier series coefficients for h
a0 <- 1/(2*pi) * integrate(h, -pi, pi)$value

for (n in 1:N) {
  an[n] <- 1/pi * integrate(function(t) h(t) * cos(n*t), -pi, pi)$value
  bn[n] <- 1/pi * integrate(function(t) h(t) * sin(n*t), -pi, pi)$value
  iOne[n] <- integrate(function(t) h(t) * cos(n*t), -pi, pi)$value
  iTwo[n] <- integrate(function(t) h(t) * sin(n*t), -pi, pi)$value
  cn[n] <- iOne[n] + i * iTwo[n]
}
a0
an[N]
bn[N]
cn[N]
#-------------------------------------------------------------------------------
#3)




