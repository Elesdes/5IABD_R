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
