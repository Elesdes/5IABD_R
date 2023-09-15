
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

# Exercice 2.5.1
# 1.2)
TF <- function(x) {
  n <- length(x)
  k <- 0:(n-1)
  omega <- 2*pi*k / n
  X <- (1/n) * x * exp(-1i*omega*k)
  return(X)
}
# 1.1)
t <- seq(-10, 10, by=0.01)
x_t <- numeric(length(t))
for (i in 1:length(t)) {
  if (abs((t[i]-pi)/(2 * pi)) <= 1/2) {
    x_t[i] <- 1
  } else {
    x_t[i] <- 0
  }
}

X_t <- TF(x_t)
plot(t, x_t, type = "l", xlab = "t", ylab = "x(t)", main="Signal X")

# Pour Y:
y_t <- 1 / sqrt(2*pi) * exp(-t^2 / 2)
Y_t <- TF(y_t)
plot(t, y_t, type = "l", xlab = "t", ylab = "y(t)", main="Signal Y")

# Pour Z f0=1:
f0 <- 1
z_t_1 <- cos(2*pi*f0*t)^2
Z_t_1 <- TF(z_t_1)
plot(t, z_t_1, type = "l", xlab = "t", ylab = "z(t) f0=1", main="Signal Z f0=1")

# Pour Z f0=2:
f0 <- 2
z_t_2 <- cos(2*pi*f0*t)^2
Z_t_2 <- TF(z_t_2)
plot(t, z_t_2, type = "l", xlab = "t", ylab = "z(t) f0=2", main="Signal Z f0=2")

# Pour Z f0=3:
f0 <- 3
z_t_3 <- cos(2*pi*f0*t)^2
Z_t_3 <- TF(z_t_3)
plot(t, z_t_3, type = "l", xlab = "t", ylab = "z(t) f0=3", main="Signal Z f0=3")

# 1.3)
plot(t, abs(X_t), main="Amplitude X")
plot(t, atan(X_t), main="Phase X")

plot(t, abs(Y_t), main="Amplitude Y")
plot(t, atan(Y_t), main="Phase Y")

plot(t, abs(Z_t_1), main="Amplitude Z f0=1")
plot(t, atan(Z_t_1), main="Phase Z f0=1")

plot(t, abs(Z_t_2), main="Amplitude Z f0=2")
plot(t, atan(Z_t_2), main="Phase Z f0=2")

plot(t, abs(Z_t_3), main="Amplitude Z f0=3")
plot(t, atan(Z_t_3), main="Phase Z f0=3")

# 2.1)
TF_r <- function(X) {
  n <- length(X)
  k <- 0:(n-1)
  omega <- 2*pi*k / n
  x <- (n) * X * exp(1i*omega*k)
  return(x)
}
# 2.2)
Exo_2 <- function(f){
  return(1/(1 + (f^2)))
}

tf1_result <- numeric(length(seq(-0.4, 0.4, by=0.01)))
iterator <- 0
for(i in seq(-0.4, 0.4, by = 0.01)){
  t_exo_2 <- i
  Integrand <- function(f) {
    return((Exo_2(f)*exp(2i * pi * f * t_exo_2)))
  }
  tf1_result[iterator] <- integral(Integrand, xmin = -Inf, xmax = Inf)
  #cat(i,":", "tf1_result(t) =", tf1_result[iterator], "\n")
  iterator <- iterator+1
}
plot(seq(-0.4, 0.4, by=0.01), abs(tf1_result), main = "Exo 2 X(f)")

# 2.3)
x_t <- TF_r(X_t)
plot(t, x_t, type = "l", xlab = "t", ylab = "x(t)", main="Signal X")
y_t <- TF_r(Y_t)
plot(t, y_t, type = "l", xlab = "t", ylab = "y(t)", main="Signal Y")
z_t_1 <- TF_r(Z_t_1)
plot(t, z_t_1, type = "l", xlab = "t", ylab = "z(t) f0=1", main="Signal Z f0=1")
z_t_2 <- TF_r(Z_t_2)
plot(t, z_t_2, type = "l", xlab = "t", ylab = "z(t) f0=2", main="Signal Z f0=2")
z_t_3 <- TF_r(Z_t_3)
plot(t, z_t_3, type = "l", xlab = "t", ylab = "z(t) f0=3", main="Signal Z f0=3")

# 3. TFTD & TFD

# Function to compute the Discrete Fourier Transform (DFT)
compute_dft <- function(x) {
  N <- length(x)
  k <- 0:(N - 1)
  Xk <- numeric(N)
  for (k_val in k) {
    Xk[k_val + 1] <- sum(x * exp(-2i * pi * k_val * (0:(N - 1)) / N))
  }
  return(Xk)
}

# Function to display a sampled signal and its amplitude spectrum
display_sampled_signal_and_amp_spectrum <- function(n, xn, Ak, fk) {
  par(mfrow=c(2, 1))
  
  # Plot sampled signal
  plot(n, xn, type="o", pch=19, xlab="n", ylab="x[n]", main="Sampled signal xn")
  
  # Plot amplitude spectrum
  plot(fk, Ak, type="o", pch=19, xlab="Frequency (Hz)", ylab="Amplitude", main="Amp spectrum")
}

# Function to compute and display the inverse Fourier Transform
compute_and_display_inverse_fourier_transform <- function(Xk) {
  x <- inverse_dft(Xk)
  display_inverse_fourier_transform(x)
}

# Function to display the inverse Fourier Transform
display_inverse_fourier_transform <- function(x) {
  x_real <- Re(x)
  plot(n, x_real, type="o", pch=19, xlab="n", ylab="x[n]", main="Signal from Inverse Fourier Transform")
}

# 3.1
fe <- 16
Te <- 1 / fe
N <- 8
  
n <- 0:(N - 1)
xn <- 2 * sin(8 * pi * n * Te) + 8 * cos(4 * pi * n * Te)
  
k <- 0:(N - 1)
Xk <- compute_dft(xn)
Ak <- abs(Xk)
fk <- k * fe / N
  
display_sampled_signal_and_amp_spectrum(n, xn, Ak, fk)

# 3.2
fe <- 16
Te <- 1 / fe
N <- 24
  
n <- 0:(N - 1)
xn <- 3 * sin(8 * pi * n * Te) + 4 * cos(6 * pi * n * Te)
  
k <- 0:(N - 1)
Xk <- compute_dft(xn)
Ak <- abs(Xk)
fk <- k * fe / N
  
display_sampled_signal_and_amp_spectrum(n, xn, Ak, fk)

# 3.3
inverse_dft <- function(Xk) {
  N <- length(Xk)
  n <- 0:(N - 1)
  x <- numeric(N)
  for (n_val in n) {
    x[n_val + 1] <- (1 / N) * sum(Xk * exp(2i * pi * n_val * (0:(N - 1)) / N))
  }
  return(x)
}

# 3.4
library(microbenchmark)

fft_algorithm <- function(x) {
  N <- length(x)
  if (N == 1) {
    return(x)
  } else {
    even_indices <- seq(2, N, by = 2)
    odd_indices <- seq(1, N, by = 2)
    even_fft <- fft_algorithm(x[even_indices])
    odd_fft <- fft_algorithm(x[odd_indices])
    twiddle_factors <- exp(-2i * pi * (0:(N / 2 - 1)) / N)
    return(c(even_fft + twiddle_factors * odd_fft, even_fft - twiddle_factors * odd_fft))
  }
}

set.seed(42)
signal <- runif(1024)

execution_time <- microbenchmark(
  FFT_result <- fft_algorithm(signal),
  times = 100
)

print(paste("Average execution time: ", summary(execution_time)$median, "ms"))
