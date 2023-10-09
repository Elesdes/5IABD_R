# Imports
library(pracma)
library(factoextra)
library(FactoMineR)
library(microbenchmark)
library(dplyr)

# 1. Série de Fourier

# 1.1
f <- function(t) { ifelse(t >= -pi & t <= pi, pi - abs(t), pi - abs(t - 2*pi*floor((t + pi)/(2*pi)))) }
g <- function(t) { ifelse(t >= -pi & t <= pi, t^2 - pi^2, (t - 2*pi*floor((t + pi)/(2*pi)))^2 - pi^2) }
h <- function(t) {
  t <- t %% (2*pi)
  ifelse(t >= 0 & t < pi, exp((-t )/pi), 0) 
}

x <- seq(-5 * pi, 5 * pi, by= 0.01)

plot(x, f(x), type = "l", main = "Graph of F", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")
plot(x, g(x), type = "l", main = "Graph of G", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")
plot(x, h(x), type = "l", main = "Graph of H", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "green")

# 1.2
# SAISIR LA VALEUR DE N
N <- as.numeric(readline("Saisissez la valeur de N :"))
#N <- 20

# CALCUL DES A0 + INITIALISATION DES AN ET BN
af0 <- 1 / (2 * pi) * integrate(f, -pi, pi)$value
afn <- rep(0, N)
bfn <- rep(0, N)

ag0 <- 1/(2*pi) * integrate(g, -pi, pi)$value
agn <- rep(0, N)
bgn <- rep(0, N)

ah0 <- 1/(2*pi) * integrate(h, -pi, pi)$value
ahn <- rep(0, N)
bhn <- rep(0, N)

# CALCUL DES AN ET BN
# F
for (n in 1:N) {
  afn[n] <- 1/pi * integrate(function(x) f(x) * cos(n*x), -pi, pi)$value 
  bfn[n] <- 1/pi * integrate(function(x) f(x) * sin(n*x), -pi, pi)$value
}
cfn <- complex(real = afn, imaginary = - bfn)
af0
afn
bfn
cfn

# G
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
# H
for (n in 1:N) {
  ahn[n] <- 1/pi * integrate(function(x) h(x) * cos(n*x), -pi, pi)$value
  bhn[n] <- 1/pi * integrate(function(x) h(x) * sin(n*x), -pi, pi)$value

}
chn <- complex(real = ahn, imaginary = -bhn)
ah0
ahn
bhn
chn

# 1.3

# SERIE REEL
sr <- function(t, an, bn, a0, N){
  rez <- 0
  for (n in 1:N) {
    rez <- rez + (an[n] * cos(n * t) + bn[n] * sin(n * t))
  }
  rez = a0 + rez
  return(rez)
}

# SERIE COMPLEXE
sc <-function(t, cn, N, a0){
  rez2 <- 0
  for (n in 1:N){
      rez2 <- rez2 + (cn[n] * exp(1i * t * n))
    }
  rez2 <- rez2 + a0
  return(rez2)
}

N <- 5

# F
srf <- sr(x, afn, bfn, af0, N)
scf <- sc(x, cfn, N, af0)

plot(x, f(x), type = "l", main = "Graph of F Reel and Complex", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "red")

lines(x, srf, type = "l", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "blue")
lines(x, scf, type = "l", xlab = "t", ylab = "f(t)", ylim = c(-20, 20), col = "green")

# G
srg <- sr(x, agn, bgn, ag0, N)
scg <- sc(x, cgn, N, ag0)

plot(x, g(x), type = "l", main = "Graph of G Reel and Complex", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "red")

lines(x, srg, type = "l", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "blue")
lines(x, scg, type = "l", xlab = "t", ylab = "g(t)", ylim = c(-20, 20), col = "green")

# H
srh <- sr(x, ahn, bhn, ah0, N)
sch <- sc(x, chn, N, ah0)

plot(x, h(x), type = "l", main = "Graph of H Reel and Complex", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "red")

lines(x, srh, type = "l", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "blue")
lines(x, sch, type = "l", xlab = "t", ylab = "h(t)", ylim = c(-20, 20), col = "green")

# 1.4
x2 <- seq(1, 14)

# AMPLITUDE
Amplin <- function(an, bn) {
  rez4 <- rep(1:14)
  for(i in 1:14) {
    rez4[i] <- sqrt(an[i] ^ 2 + bn[i] ^ 2) / 2
  }
  return(rez4)
}

# PHASE
phaseN <- function(N, t, func){
  rez3 <- 0
    for(n in 1:N){
      rez3 <- rez3 + (exp((2 * n + 1) * 1i * t) / (2 * n + 1))
    }
    rez3 <- ((-2i * func) / pi) * rez3
  return(rez3)
}

Amplifn <- Amplin(afn, bfn) 
Amplign <- Amplin(agn, bgn)
Amplihn <- Amplin(ahn, bhn)

t <- 1
phaseNf <- phaseN(14, t, f(t))
phaseNg <- phaseN(14, t, g(t))
phaseNh <- phaseN(14, t, h(t))

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

# 2. TF

# 2.1
TF <- function(x) {
  n <- length(x)
  k <- 0:(n - 1)
  omega <- 2 * pi * k / n
  X <- (1 / n) * x * exp(-1i * omega * k)
  return(X)
}

t <- seq(-10, 10, by=0.01)
x_t <- numeric(length(t))
for (i in 1:length(t)) {
  if (abs((t[i] - pi)/(2 * pi)) <= 1 / 2) {
    x_t[i] <- 1
  } else {
    x_t[i] <- 0
  }
}

X_t <- TF(x_t)
plot(t, x_t, type = "l", xlab = "t", ylab = "x(t)", main="Signal X")

# Pour Y:
y_t <- 1 / sqrt(2 * pi) * exp(-t ^ 2 / 2)
Y_t <- TF(y_t)
plot(t, y_t, type = "l", xlab = "t", ylab = "y(t)", main="Signal Y")

# Pour Z f0=1:
f0 <- 1
z_t_1 <- cos(2 * pi * f0 * t) ^ 2
Z_t_1 <- TF(z_t_1)
plot(t, z_t_1, type = "l", xlab = "t", ylab = "z(t) f0=1", main="Signal Z f0=1")

# Pour Z f0=2:
f0 <- 2
z_t_2 <- cos(2 * pi * f0 * t) ^ 2
Z_t_2 <- TF(z_t_2)
plot(t, z_t_2, type = "l", xlab = "t", ylab = "z(t) f0=2", main="Signal Z f0=2")

# Pour Z f0=3:
f0 <- 3
z_t_3 <- cos(2 * pi * f0 * t)^2
Z_t_3 <- TF(z_t_3)
plot(t, z_t_3, type = "l", xlab = "t", ylab = "z(t) f0=3", main="Signal Z f0=3")

# Display amp and phase spectrum
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

# Reverse TF
TF_r <- function(X) {
  n <- length(X)
  k <- 0:(n - 1)
  omega <- 2 * pi * k / n
  x <- (n) * X * exp(1i * omega * k)
  return(x)
}

# 2.2
Exo_2 <- function(f){
  return(1 / (1 + (f ^ 2)))
}

tf1_result <- numeric(length(seq(-0.4, 0.4, by=0.01)))
iterator <- 0
for(i in seq(-0.4, 0.4, by = 0.01)){
    t_exo_2 <- i
    Integrand <- function(f) {
      return((Exo_2(f) * exp(2i * pi * f * t_exo_2)))
    }
    tf1_result[iterator] <- integral(Integrand, xmin = -Inf, xmax = Inf)
    iterator <- iterator + 1
  }
plot(seq(-0.4, 0.4, by=0.01), abs(tf1_result), main = "Exo 2 X(f)")

# 2.3
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

print(paste("The average execution time of the FFT algorithm is", summary(execution_time)$median, "ms."))

# Part 4 Méthode du Chi deux
# Le nombre total de client est de 200. Prenons cette base.

theorical <- c(0.1 * 200, 0.1 * 200, 0.15 * 200, 0.2 * 200, 0.3 * 200, 0.15 * 200)
obs <- c(30, 14, 34, 45, 57, 20)

khi_deux <- numeric(length(6))
for (i in 0:6) {
  khi_deux[i] <- (obs[i]-theorical[i])^2/theorical[i]
}
khi_deux_total <- sum(khi_deux)
cat("Khi deux: ", khi_deux, "\n")
cat("Khi deux total: ", khi_deux_total)
k <- 5
alpha <- 0.05
Q <- 11.07

# Malheureusement, Khi deux total est supérieur à la valeur du seuil (11.4>11.07)
# Nous ne devrions pas ouvrir le restaurant

# Part 5 ACP
getwd()

# A)
donnees <- read.table("data/input/decathlon.dat", header=TRUE, sep=" ")
sous_donnees <- donnees[, 1:(ncol(donnees)-1)]
sous_donnees <- apply(sous_donnees, 2, as.numeric)

matrice_correlation <- cor(sous_donnees)
print(matrice_correlation)
write.csv2(matrice_correlation, file = "data/output/matrice_correlation.csv")

# Il est assez dur d'interpréter les résultats.
# Il semblerait que le nombre de Points est fortement corrélé à l'épreuve du saut en longueur,
# du lancer de poids et du saut en hauteur.
# On remarque que le lancer de disque est souvent corrélé au lancer de poids.
# Les trois courses (100, 110 et 400m) sont bien corrélées.
# Les épreuves de saut en longueur sont fortement opposées au rang et aux trois courses.
# Les points sont fortement opposés aux trois courses.
# Enfin, le reste de ce qui n'a pas été cité est très souvent peu corrélé.

# Plus les variables ont un ratio proche de 1 et plus elles sont corrélées entre elles.
# Cela signifie que si l'une augmente, l'autre à tendance à augmenter aussi.
# En revanche, plus elles sont proches de -1 et plus l'une augmente, plus l'autre diminue.
# Cela reste une corrélation bien que négative. Cela s'explique simplement par un calcul de multiplication:
# Si A est positif et que le signe se rapproche de -1, alors B sera de plus en plus petit.

# B)
data <- read.table("data/input/decathlon.dat", header=TRUE, sep=" ")
sub_data <- data[, 1:(ncol(data)-3)]
sub_data <- apply(sub_data, 2, as.numeric)
mat_cor <- cor(sub_data)
# Analyse en composantes principales
res.pca <- prcomp(mat_cor, scale = TRUE)
fviz_eig(res.pca)

eigenvalues <- get_eigenvalue(res.pca)
print("Valeurs propres")
print(eigenvalues)
eigenvalues <- sort(eigenvalues$eigenvalues, decreasing=TRUE)
# On constate qu'au bout de 3 dimensions, on obtient une variance cumulée à plus de 87%.
# Fixons ainsi notre seuil à 3 dimensions.
# En effet, en regardant la courbe, le but est de trouver le critère de Kaiser appelé
# aussi la règle du coude. Cette règle peut nous emmener jusqu'à la 4ème dimension si nous
# le souhaitons mais par soucis de praticité, nous allons nous arrêter à 3.
# Nous avons donc une inertie à 87%.
pca_result <- PCA(mat_cor, graph=FALSE, ncp=3)
C1 <- pca_result$ind$coord[,1]
C2 <- pca_result$ind$coord[,2]
C3 <- pca_result$ind$coord[,3]

print(C1)
print(C2)
print(C3)
# Il est clair et net que la C1 correspond principalement aux courses et le saut en longueur,
# le lancer de poids et le lancer de disque.
# C2 correspond presque uniquement au saut à la perche et un peu la course sur 1500m.
# C3 va plus concerner la corrélation entre le saut en longueur/à la perche, la course sur 1500m et le lancer de disque/javelot

correlation_C1 <- cor(mat_cor, C1)
correlation_C2 <- cor(mat_cor, C2)
correlation_C3 <- cor(mat_cor, C3)

cat("Tableau de correlation C1: \n", round(correlation_C1,2), "\n")
cat("Tableau de correlation C2: \n", round(correlation_C2,2), "\n")
cat("Tableau de correlation C3: \n", round(correlation_C3,2), "\n")

circle_C1_C2 <- fviz_pca_var(pca_result, axes = c(1, 2), col.var = "black", title = "Cercle de corrélation (C1, C2)")
circle_C2_C3 <- fviz_pca_var(pca_result, axes = c(2, 3), col.var = "black", title = "Cercle de corrélation (C2, C3)")

print(circle_C1_C2)
print(circle_C2_C3)

correlations_table <- data.frame(Variable = colnames(mat_cor), C1 = correlation_C1, C2 = correlation_C2, C3 = correlation_C3)
print(correlations_table)

# Il est communément admis que 0.5 est un seuil convenable comme 0.3. Cependant puisque nous avons des
# variables fortement corrélées, nous allons prendre 0.5.

seuil <- 0.5
for (i in 1:3) {
  cercle_correlation <- pca_result$var$coord
  variables_influentes <- rownames(cercle_correlation)[abs(cercle_correlation[, i]) > seuil]
  cat("Variables influentes pour C", i, ":", paste(variables_influentes, collapse = ", "), "\n")
}

# Il est très dur de savoir laquelle des variables est la plus influentes pour C1 entre C100, C110 et C400
# Selon les chiffres, C110 serait la composante principale.
# Pour C2, il est clair et net qu'il s'agit bien plus du saut à la perche que le lancer de javelot et C1500
# Pour C3, il semblerait qu'il s'agisse bien plus du C1500 qui joue.

# L'effet de taille est déjà présent puisque l'on applique un effet de centré-réduit pour mettre à la même échelle
# les valeurs.




#PARTIE 6 - ACP
df <- read.csv("data/input/donnees_sympathiques.csv", header = TRUE)

# Exclure la colonne TOTAL, la colonne CAT, et la dernière ligne qui contient les totaux
characteristics <- df[-nrow(df), -c(1, ncol(df))]

# Somme pour chaque caractéristique
sum_characteristics <- colSums(characteristics)

# Total des votes
total_votes <- sum(sum_characteristics)

# Pourcentage pour chaque caractéristique
percentage_characteristics <- (sum_characteristics / total_votes) * 100
print("Pourcentage pour chaque caractéristique:")
print(percentage_characteristics)

# Somme pour chaque catégorie professionnelle
sum_categories <- rowSums(characteristics)

# Pourcentage pour chaque catégorie
percentage_categories <- (sum_categories / total_votes) * 100
print("Pourcentage pour chaque catégorie professionnelle:")
names(percentage_categories) <- df$CAT[-nrow(df)] 
print(percentage_categories)

#Question 2
df <- read.csv("data/input/donnees_sympathiques.csv", header = TRUE)
df

# TOTAL des personnes sondées
total_people <- df[nrow(df), "TOTAL"] / 3

print(paste("Total de personnes sondées :", total_people))

# Proportion des employés (EMPL) pour qui être honnête rend sympathique
employees_honest <- df[df$CAT == "EMPL", "HONN"]
total_employees <- df[df$CAT == "EMPL", "TOTAL"]
percentage_employees_honest <- (employees_honest / total_employees) * 100
print(paste("Proportion des employés pour qui être honnête rend sympathique :", round(percentage_employees_honest, 2), "%"))


# Proportion d'employés (EMPL) parmi les gens qui pensent qu'être honnête rend sympathique
total_honest <- df[nrow(df), "HONN"]
percentage_employees_among_honest <- (employees_honest / total_honest) * 100
print(paste("Proportion d'employés parmi les gens qui pensent qu'être honnête rend sympathique :", round(percentage_employees_among_honest), "%"))

#-----------------------------------------------------------------------------------------------------------------------------------

# Création du tableau de contingence (exemple fictif)
# Exemple de tableau de contingence correct
tableau_contingence <- matrix(c(20, 9, 9, 27, 10, 16, 20, 4, 8,
                                42, 10, 22, 51, 18, 28, 38, 12, 22,
                                11, 2, 5, 14, 8, 7, 5, 8, 6,
                                8, 9, 12, 23, 14, 16, 14, 12, 12,
                                19, 10, 16, 52, 32, 25, 22, 25, 30,
                                10, 5, 12, 23, 20, 13, 11, 13, 10,
                                2, 8, 7, 6, 15, 6, 6, 9, 4,
                                8, 42, 23, 24, 46, 22, 22, 34, 16),
                              nrow = 8, ncol = 9, byrow = TRUE)

# Noms des variables (colonnes)
colnames(tableau_contingence) <- c("SERI", "GENE", "GAI", "HONN", "INTL", "SERV", "COUR", "COMP", "DISC")

# Nom des lignes (catégories)
rownames(tableau_contingence) <- c("PAYS", "OUVR", "VEND", "COMM", "EMPL", "TECH", "UNIV", "LIBE")

# Effectuer l'Analyse Factorielle des Correspondances
afc_result <- CA(tableau_contingence)

# Obtenir les valeurs propres
valeurs_propres <- afc_result$eig
print(valeurs_propres)
