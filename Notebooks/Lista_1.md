MI628 - Lista 1
================

# Pacotes e Funções

``` r
# PACOTES
library(Matching)
```

    ## Carregando pacotes exigidos: MASS

    ## ## 
    ## ##  Matching (Version 4.10-14, Build Date: 2023-09-13)
    ## ##  See https://www.jsekhon.com for additional documentation.
    ## ##  Please cite software as:
    ## ##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching
    ## ##   Software with Automated Balance Optimization: The Matching package for R.''
    ## ##   Journal of Statistical Software, 42(7): 1-52. 
    ## ##

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ dplyr::select() masks MASS::select()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# FUNCOES
# Função que ajusta modelos para todas as combinações de preditoras
# possíveis
ajuste.subconjuntos <- function(df) {
  
  colunas <- colnames(lalonde)[colnames(lalonde) != c("re78", "treat")]
  
  combin <- expand.grid(v1 = c(T,F), v2 = c(T,F), v3 = c(T,F), v4 = c(T,F),
                      v5 = c(T,F), v6 = c(T,F), v7 = c(T,F), v8 = c(T,F),
                      v9 = c(T,F), v10 = c(T,F))
  colnames(combin) <- colunas
  l <- unlist(apply(combin, 1, list), recursive = FALSE)
  subconjuntos <- lapply(l, function(x) names(x)[x])
  
  ajustes <- list()
  
  for (ajuste in seq_len(length(subconjuntos))) {
    
    dados <- df %>% 
      select(re78, treat, subconjuntos[[ajuste]]) %>% 
      na.omit()
    
    if (length(dados) == 2) {
      
      fit <- formula(paste0("re78 ~ treat"))
      
    } else {
      
      fit <- formula(paste0("re78 ~ treat + ", 
                          paste(subconjuntos[[ajuste]], collapse = " + ")))
      
    }
    
    modelo <- lm(formula = fit, data = dados)
    
    ajustes[[ajuste]] <- modelo
    
  }

  return(ajustes)
  
}

# Função para extrair todos os betas da lista de ajustes
extraindo.betas <- function(ajustes) {
  
  betas <- c()
  
  for (ajuste in seq_len(length(ajustes))) {
    
    beta <- ajustes[[ajuste]]$coefficients[2]
    
    betas <- c(betas, beta)
    
  }
  
  return(betas)
  
}

# Função para extrair todos os p-valores da lista de ajustes
extraindo.pvalor <- function(ajustes) {
  
  pvalores <- c()
  
  for (ajuste in seq_len(length(ajustes))) {
    
    pvalor <- summary(ajustes[[ajuste]])$coefficients[2, 3]
    
    pvalores <- c(pvalores, pvalor)
    
  }
  
  return(pvalores)
  
}
```

# Questão 3

``` r
# Carregando base
data(lalonde)

# Transformando variável treat em dummy
lalonde <- lalonde %>% 
  mutate(treat = factor(treat))
```

``` r
# Ajustando 1024 regressões
ajustes <- ajuste.subconjuntos(lalonde)

# Extraindo os estimadores para treat
betas <- extraindo.betas(ajustes)

# Extraindo p-valores dos estimadores
pvalores <- extraindo.pvalor(ajustes)

# Salvando resultados
resultados <- tibble(betas, pvalores) %>% 
  mutate(pos_sig = ifelse(betas > 0 & pvalores > 2.807034, 1, 0),
         neg_sig = ifelse(betas < 0 & pvalores > 2.807034, 1, 0),
         sig = ifelse(pvalores <= 2.807034, 1, 0))
```

1)  Quantidade de vezes em que o tratamento foi positivo e
    significativo: 106
2)  Quantidade de vezes em que o tratamento foi negativo e
    significativo: 0
3)  Quantidade de vezes em que o tratamento não foi significativo: 918

``` r
library(car)
```

    ## Carregando pacotes exigidos: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

    ## The following object is masked from 'package:purrr':
    ## 
    ##     some

``` r
library(purrr)

se_EHW <- map(ajustes, function(x) sqrt(hccm(x, type = "hc0")[2, 2]))
Z <- map2(.x = ajustes, .y = se_EHW, .f = function(x, y) x$coefficients[2]/y)


n_signficativos <- as.numeric(unlist(map(Z, function(x) ifelse(x > 2.807034, 1, 0))))
sum(n_signficativos)
```

    ## [1] 1

# Questão 4a

``` r
E_rho_pos <- -0.5 - 2 * ((-dnorm(0.5))/(1 - pnorm(0.5, lower.tail = T)))

E_rho_neg <- - 1 * (dnorm(0.5)/(pnorm(0.5)))

E_diff <- E_rho_pos - E_rho_neg

E_diff
```

    ## [1] 2.291316

# Questão 4b

``` r
Y0 <- rnorm(1)
tau <- -0.5 + Y0
Y1 <- Y0 + tau
Z <- ifelse(tau < 0, 0, 1)
Y <- Z * Y1 + (1-Z) * Y0

Y
```

    ## [1] -0.5604756

``` r
Y0 <- rnorm(1000000)
tau <- -0.5 + Y0
Y1 <- Y0 + tau
Z <- ifelse(tau < 0, 0, 1)
Y <- Z * Y1 + (1-Z) * Y0

mean(Y[Z==1])
```

    ## [1] 1.780949

``` r
mean(Y[Z==0])
```

    ## [1] -0.5096613

``` r
mean(Y[Z==1]) - mean(Y[Z==0])
```

    ## [1] 2.29061
