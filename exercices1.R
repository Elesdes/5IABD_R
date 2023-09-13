#1) Analyse de Fourier




f <- function(t) { ifelse(t >= -pi & t <= pi, pi - abs(t), pi - abs(t - 2*pi*floor((t + pi)/(2*pi)))) }


g <- function(t) { ifelse(t >= -pi & t <= pi, t^2 - pi^2, (t - 2*pi*floor((t + pi)/(2*pi)))^2 - pi^2) }



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

#N <- as.numeric(readline("Saisissez la valeur de N :"))
N <- 20



af0 <- 1/(2*pi) * integrate(f, -pi, pi)$value
afn <- rep(0, N)
bfn <- rep(0, N)
ifOne <- rep(0, N)
ifTwo <- rep(0, N)
cfn <- rep(0, N)


ag0 <- 1/(2*pi) * integrate(f, -pi, pi)$value
agn <- rep(0, N)
bgn <- rep(0, N)
igOne <- rep(0, N)
igTwo <- rep(0, N)
cgn <- rep(0, N)


ah0 <- 1/(2*pi) * integrate(f, -pi, pi)$value
ahn <- rep(0, N)
bhn <- rep(0, N)
ihOne <- rep(0, N)
ihTwo <- rep(0, N)
chn <- rep(0, N)



for (n in 1:N) {
  afn[n] <- 1/pi * integrate(function(x) f(x) * cos(n*x), -pi, pi)$value
  bfn[n] <- 1/pi * integrate(function(x) f(x) * sin(n*x), -pi, pi)$value
  ifOne[n] <- integrate(function(x) f(x) * cos(n*x), -pi, pi)$value
  ifTwo[n] <- integrate(function(x) f(x) * sin(n*x), -pi, pi)$value
  cfn[n] <- ifOne[n] + (-1i) * ifTwo[n]
}
af0
afn
bfn
cfn
# Calculate the Fourier series coefficients for g
for (n in 1:N) {
  agn[n] <- 1/pi * integrate(function(x) g(x) * cos(n*x), -pi, pi)$value
  bgn[n] <- 1/pi * integrate(function(x) g(x) * sin(n*x), -pi, pi)$value
  igOne[n] <- integrate(function(x) g(x) * cos(n*x), -pi, pi)$value
  igTwo[n] <- integrate(function(x) g(x) * sin(n*x), -pi, pi)$value
  cgn[n] <- igOne[n] + 1i * igTwo[n]
}
ag0
agn[N]
bgn[N]
cgn[N]

# Calculate the Fourier series coefficients for h


for (n in 1:N) {
  ahn[n] <- 1/pi * integrate(function(x) h(x) * cos(n*x), -pi, pi)$value
  bhn[n] <- 1/pi * integrate(function(x) h(x) * sin(n*x), -pi, pi)$value
  ihOne[n] <- integrate(function(x) h(x) * cos(n*x), -pi, pi)$value
  ihTwo[n] <- integrate(function(x) h(x) * sin(n*x), -pi, pi)$value
  chn[n] <- ihOne[n] + 1i * ihTwo[n]
}
ah0
ahn
bhn
chn
#-------------------------------------------------------------------------------
#3)
N <- 1

sr <- function(t, an, bn){
  rez <- 0
  for (n in 1:N){
    rez <- rez + (an[n] * cos(n*t) + bn[n] * sin(n * t))
  }
  rez = af0 + rez
  return(rez)
}

sc <-function(t, cn){
  rez2 <- 0
  for (n in 1:N){
    rez2 <- rez2 + (cn[n] * exp(1i*t))
  }
  
  return(rez2)
}

srf <- sr(x, afn, bfn)
scf <- sc(x, cfn)
plot(x, f(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")

lines(x, srf, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "blue")
lines(x, scf, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "black")




srg <- sr(x, agn, bgn)
scg <- sc(x, cgn)

plot(x, g(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "red")

lines(x, srg, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")
lines(x, scg, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "black")




srh <- sr(x, ahn, bhn)
srh
sch <- sc(x, chn)

plot(x, h(x), type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "red")

lines(x, srh, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "blue")
lines(x, sch, type = "l", main = "Graph of Periodic Functions", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "black")




















