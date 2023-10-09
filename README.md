# Mathématiques pour le Big Data - Projet Semaine Thématique

## Setup Project

Packages à installer : `install.packages(c("pracma", "microbenchmark", "nycflights13", "factoextra", "dplyr", "magrittr", "FactoMineR"))`

## Projet : Analyse de Fourier, ACP et AFC

### Consignes

- Résoudre les différents exercices à l'aide du logiciel R sous R-studio.
- Écrire les scripts pour chaque traitement.
- Présenter les résultats dans un PDF avec tous les tracés et tableaux de calcul.
- Commenter et interpréter les résultats lors de la soutenance.

## Série de Fourier (SF)

### Formulation générale

Une fonction réelle f (ou signal) de période T, C1 par morceaux sur tout intervalle [α, α + T], α ∈ ℝ quelconque, se développe en une série de Fourier de la forme :
```
f(t) = a0 + ∑ [an cos(nωt) + bn sin(nωt)]
       n=1
```
avec ω la pulsation du signal.

**NB :** Aux points de discontinuité, f(t) est remplacée par 1/2 (f(t+) + f(t-)).

### Vocabulaire

1. Une fonction périodique se répète dans le temps. Sa période est mesurée en secondes (T) et la fréquence en Hertz (F = 1/T).
2. ω = 2π/T : la pulsation du signal en rad/s.
3. f est C1 par morceaux sur [a, b] : f est continue, dérivable, de dérivée continue sur [a, b] sauf en un nombre fini de points.

### Forme complexe des SF

La série de Fourier peut également être exprimée en forme complexe :
```
f(t) = ∑ cne^(inωt), avec cn = (1/T) ∫[α, α+T] f(t)e^(-inωt)dt
```

### Harmonique, spectre de fréquence

1. Pour la série de Fourier réelle : Le spectre d'amplitude est représenté par un diagramme montrant l'amplitude des composantes en fonction de leur rang (n).

2. Pour la série de Fourier complexe : Le spectre en fréquence est symétrique par rapport à l'axe des fréquences. La composante n = 1 est le fondamental, et les autres composantes sont les harmoniques.

### Exercice : études de cas

Dans cet exercice, vous étudierez trois signaux périodiques et calculerez leur série de Fourier, représenterez les sommes des séries réelles et complexes tronquées, ainsi que les spectres d'amplitude et de phase.

## Transformée de Fourier (TF)

### Intuitivement

La transformée de Fourier permet de passer du domaine temporel/spatial au domaine fréquentiel/spectral. Elle permet de décomposer un signal en ses différentes composantes fréquentielles.

### Formules

La transformée de Fourier d'un signal x(t) est définie comme suit :
```
X(f) = ∫ x(t)e^(-2πift)dt
```

### Spectre de fréquence

On peut obtenir le spectre d'amplitude et le spectre de phase à partir de la transformée de Fourier complexe X(f).

## Transformée de Fourier Discrète (TFD)

La TFD est utilisée pour analyser les signaux numériques discrets. Elle permet de calculer une représentation spectrale discrète d'un signal échantillonné sur une fenêtre de temps finie.

### Formules

La TFD d'un signal de N échantillons (x(0), x(1), ..., x(N-1)) est définie comme suit :
```
X(k) = Σ x(n)e^(-2πikn/N), 0 <= k < N
```

### Exercices

Dans cette section, vous résoudrez des exercices liés à la TFD, à la TFD inverse et à l'algorithme FFT.

## Méthode du Chi-deux

La méthode du Chi-deux est utilisée pour déterminer la nature d'une répartition. Elle compare les valeurs observées et les valeurs théoriques attendues pour évaluer le désaccord entre les deux.

### Exemple

Dans cet exemple, nous utilisons la méthode du Chi-deux pour déterminer si un dé est pipé ou non en comparant les résultats observés et les résultats théoriques attendus.

## Analyse en Composantes Principales (ACP)

L'ACP est une technique de réduction de dimensionnalité qui permet d'explorer la structure d'un ensemble de données en identifiant les directions principales de variation.

### Étapes de l'ACP

1. Standardisation des données.
2. Calcul de la matrice de covariance.
3. Calcul des vecteurs propres et des valeurs propres.
4. Sélection des composantes principales.
5. Projection des données sur les composantes principales.

### Exercices

Dans cette section, vous appliquerez l'ACP à un jeu de données de dimension supérieure et interpréterez les résultats.

## Analyse Factorielle des Correspondances (AFC)

L'AFC est une technique d'analyse des données qui permet d'explorer les relations entre deux ensembles de catégories. Elle est utilisée pour analyser des données catégorielles.

### Étapes de l'AFC

1. Construction du tableau de contingence.
2. Calcul des fréquences marginales.
3. Calcul des fréquences conditionnelles.
4. Calcul des statistiques de l'AFC.
5. Représentation des résultats.

### Exercices

Dans cette section, vous appliquerez l'AFC à un jeu de données catégorielles et interpréterez les résultats.

---

N'oubliez pas de commenter vos scripts de manière appropriée et de présenter vos résultats de manière claire dans le PDF de la soutenance. Bonne chance avec votre projet !