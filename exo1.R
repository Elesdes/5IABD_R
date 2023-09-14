#EXERCICE 1 - SERIE DE FOURIER 
#1) Analyse de Fourier



#FONCTION F
f <- function(t) { ifelse(t >= -pi & t <= pi, pi - abs(t), pi - abs(t - 2*pi*floor((t + pi)/(2*pi)))) }

#FONCTION G
g <- function(t) { ifelse(t >= -pi & t <= pi, t^2 - pi^2, (t - 2*pi*floor((t + pi)/(2*pi)))^2 - pi^2) }

#FONCTION H
h <- function(t) {
  t <- t %% (2*pi)
  ifelse(t >= 0 & t < pi, exp((-t )/pi), 0) 
}
x <- seq(-5 * pi, 5 * pi, by= 0.01)

plot(x, f(x), type = "l", main = "Graph of F", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")


plot(x, g(x), type = "l", main = "Graph of G", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")

plot(x, h(x), type = "l", main = "Graph of H", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "green")
#-----------------------------------------------------------------------------------------------------------------------------------
#2)


#SAISIR LA VALEUR DE N
#N <- as.numeric(readline("Saisissez la valeur de N :"))
N <- 20


#CALCUL DES A0 + INITIALISATION DES AN ET BN
af0 <- 1/(2*pi) * integrate(f, -pi, pi)$value
afn <- rep(0, N)
bfn <- rep(0, N)


ag0 <- 1/(2*pi) * integrate(g, -pi, pi)$value
agn <- rep(0, N)
bgn <- rep(0, N)


ah0 <- 1/(2*pi) * integrate(h, -pi, pi)$value
ahn <- rep(0, N)
bhn <- rep(0, N)



#CALCUL DES AN ET BN
#F
for (n in 1:N) {
  afn[n] <- 1/pi * integrate(function(x) f(x) * cos(n*x), -pi, pi)$value 
  bfn[n] <- 1/pi * integrate(function(x) f(x) * sin(n*x), -pi, pi)$value
}
cfn <- complex(real = afn, imaginary = - bfn)
af0
afn
bfn
cfn

#G
for (n in 1:N) {
  agn[n] <- 1/pi * integrate(function(x) g(x) * cos(n*x), -pi, pi)$value
  bgn[n] <- 1/pi * integrate(function(x) g(x) * sin(n*x), -pi, pi)$value
}
cgn <- complex(real = agn, imaginary = - bgn)
ag0
agn
bgn
cgn

# Calculate the Fourier series coefficients for h

#H
for (n in 1:N) {
  ahn[n] <- 1/pi * integrate(function(x) h(x) * cos(n*x), -pi, pi)$value
  bhn[n] <- 1/pi * integrate(function(x) h(x) * sin(n*x), -pi, pi)$value

}
chn <- complex(real = ahn, imaginary = -bhn)
ah0
ahn
bhn
chn
#-------------------------------------------------------------------------------
#3)

#SERIE REEL
sr <- function(t, an, bn, a0, N){
  rez <- 0
  for (n in 1:N){
    rez <- rez + (an[n] * cos(n*t) + bn[n] * sin(n * t))
  }
  rez = a0 + rez
  return(rez)
}
#SERIE COMPLEXE
sc <-function(t, cn, N, a0){
  rez2 <- 0
  for (n in 1:N){
      rez2 <- rez2 + (cn[n] * exp(1i*t*n))
    }
  rez2 <- rez2 + a0
  return(rez2)
}

N <- 5

#F
srf <- sr(x, afn, bfn, af0, N)
scf <- sc(x, cfn, N, af0)
plot(x, f(x), type = "l", main = "Graph of F Reel and Complex", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")

lines(x, srf, type = "l", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "blue")
lines(x, scf, type = "l", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "green")



#G
srg <- sr(x, agn, bgn, ag0, N)
scg <- sc(x, cgn, N, ag0)

plot(x, g(x), type = "l", main = "Graph of G Reel and Complex", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "red")

lines(x, srg, type = "l", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")
lines(x, scg, type = "l", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "green")



#H
srh <- sr(x, ahn, bhn, ah0, N)
sch <- sc(x, chn, N, ah0)

plot(x, h(x), type = "l", main = "Graph of H Reel and Complex", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "red")

lines(x, srh, type = "l", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "blue")
lines(x, sch, type = "l", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "green")



#---------------------------------------------------------------------------------------------------------------------------
#4)
x2 <- seq(1, 14)
#Amplitude 
Amplin <- function(an, bn){
  rez4 <- rep(1:14)
  for(i in 1:14){
    rez4[i] <- sqrt(an[i]^2 + bn[i]^2) / 2
  }
  return(rez4)
}



#PHASE


phaseN <- function(N, t, functionName){
  rez3 <- 0
  func = switch(functionName, "f"= f(t), "g" = g(t), "h" = h(t))  
    for(n in 1:N){
      rez3 <- rez3 + (exp((2*n + 1)*1i*t) / (2*n + 1))
    }
    rez3 <- ((-2i*func) / pi) * rez3
  return(rez3)
}



Amplifn <- Amplin(afn, bfn) 
Amplign <- Amplin(agn, bgn)
Amplihn <- Amplin(ahn, bhn)

phaseNf <- phaseN(14, 1, "f")
phaseNg <- phaseN(14, 1, "g")
phaseNh <- phaseN(14, 1, "h")

length(Amplifn)
x3 <- seq(1:14)

plot(x3, Amplifn, type = "l")
plot(x3, Amplign, type = "l")
plot(x3, Amplihn, type = "l")

phaseNfDisplay <- rep(phaseNf, 7)
phaseNgDisplay <- rep(phaseNg, 7)
phaseNhDisplay <- rep(phaseNh, 7)

plot(1:7, phaseNfDisplay)
plot(1:7, phaseNgDisplay)
plot(1:7, phaseNhDisplay)











