MI628 - Lista 3
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

# Questão 2

## Importando dados

``` r
# Importando dados
df_multicenter <- read_csv("https://dataverse.harvard.edu/api/access/datafile/7440246")
```

    ## Rows: 29 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (10): center, n0, mean0, sd0, n1, mean1, sd1, n5, mean5, sd5
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
glimpse(df_multicenter)
```

    ## Rows: 29
    ## Columns: 10
    ## $ center <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, …
    ## $ n0     <dbl> 7, 11, 6, 10, 10, 6, 5, 12, 8, 9, 15, 8, 12, 9, 6, 14, 13, 15, …
    ## $ mean0  <dbl> 0.43, 0.10, 2.58, -2.30, 2.08, 1.13, 1.20, -1.21, 1.13, -0.11, …
    ## $ sd0    <dbl> 4.58, 4.21, 4.80, 3.86, 6.46, 3.24, 7.85, 2.66, 5.28, 3.62, 6.1…
    ## $ n1     <dbl> 7, 11, 6, 10, 10, 5, 6, 13, 8, 10, 14, 9, 12, 9, 7, 14, 13, 13,…
    ## $ mean1  <dbl> -5.43, -2.59, -3.94, -1.23, -6.70, 3.40, -3.67, 0.18, -2.19, -2…
    ## $ sd1    <dbl> 5.53, 3.95, 4.25, 5.17, 7.45, 8.17, 4.89, 3.81, 5.17, 5.35, 5.3…
    ## $ n5     <dbl> 8, 12, 7, 10, 10, 5, 5, 12, 9, 10, 15, 9, 11, 7, 6, 12, 13, 14,…
    ## $ mean5  <dbl> -2.63, -2.21, 1.29, -1.40, -5.13, -1.59, -1.40, -4.08, -1.96, 0…
    ## $ sd5    <dbl> 3.38, 4.14, 7.39, 2.27, 3.91, 3.19, 2.61, 6.32, 5.84, 3.53, 4.2…

``` r
# Tamanho controle
n_controle <- sum(df_multicenter$n0)

# Tamanho tratamento 1mg
n_1mg <- sum(df_multicenter$n1)

# Tamanho tratamento 5mg
n_5mg <- sum(df_multicenter$n5)
```

## Finasterida 1mg vs Controle

``` r
# Total com tratamento 1mg
total_1mg <- n_controle + n_1mg

# Tratando base tratamento 1mg
df_multicenter_1mg <- df_multicenter %>% 
  select(colnames(df_multicenter)[!str_detect(colnames(df_multicenter), "5$")]) %>% 
  mutate(pi_k = (n0+n1)/total_1mg,
         tau_k = mean1 - mean0,
         var_k = sd1^2 / n1 + sd0^2 / n0)

# Efeito causal médio
tau_s_1mg <- sum(df_multicenter_1mg$pi_k * df_multicenter_1mg$tau_k)

# Variância
v_s_1mg <- sum(df_multicenter_1mg$pi_k^2 * df_multicenter_1mg$var_k)

# IC 95% (conservador)
lim_inf_1mg <- tau_s_1mg - 1.96*sqrt(v_s_1mg)
lim_sup_1mg <- tau_s_1mg + 1.96*sqrt(v_s_1mg)

est_s_1mg <- tau_s_1mg/sqrt(v_s_1mg)

# Teste bilateral
list("Neyman IC" = c(lim_inf_1mg, lim_sup_1mg),
     "Neyman p-valor" = 2 - 2 * pnorm(abs(est_s_1mg)))
```

    ## $`Neyman IC`
    ## [1] -1.468874  0.154194
    ## 
    ## $`Neyman p-valor`
    ## [1] 0.1123782

## Finasterida 5mg vs Controle

``` r
# Total com tratamento 5mg
total_5mg <- n_controle + n_5mg

# Tratando base tratamento 5mg
df_multicenter_5mg <- df_multicenter %>% 
  select(colnames(df_multicenter)[!str_detect(colnames(df_multicenter), "1$")]) %>% 
  mutate(pi_k = (n0+n5)/total_5mg,
         tau_k = mean5 - mean0,
         var_k = sd5^2 / n5 + sd0^2 / n0)

# Efeito causal médio
tau_s_5mg <- sum(df_multicenter_5mg$pi_k * df_multicenter_5mg$tau_k)

# Variância
v_s_5mg <- sum(df_multicenter_5mg$pi_k^2 * df_multicenter_5mg$var_k)

# IC 95% (conservador)
lim_inf_5mg <- tau_s_5mg - 1.96*sqrt(v_s_5mg)
lim_sup_5mg <- tau_s_5mg + 1.96*sqrt(v_s_5mg)

# Teste Bilateral
est_s_5mg <- tau_s_5mg/sqrt(v_s_5mg)

list("Neyman IC" = c(lim_inf_5mg, lim_sup_5mg),
     "Neyman p-valor" = 2 - 2 * pnorm(abs(est_s_5mg)))
```

    ## $`Neyman IC`
    ## [1] -2.4259519 -0.8555447
    ## 
    ## $`Neyman p-valor`
    ## [1] 4.211066e-05

# Questão 3

## Importando e tratando dados

``` r
data("lalonde")

# Criando variável race
lalonde <- lalonde %>% 
  mutate(race = ifelse(black == 1, "black",
                       ifelse(hisp == 1, "hisp", "other")))

# Vetor tratamento
z <- lalonde$treat

# Covariável race
x1 <- lalonde$race

# Resultados 
y <- lalonde$re78

n <- nrow(lalonde) # Tamanho total
n1 <- sum(z) # Tamanho tratamento
n0 <- n - n1 # Tamanho controle
tau <- mean(y[z == 1]) - mean(y[z == 0]) # Efeito causal médio
s2 <- var(y) # Variância
```

## FRT

``` r
# FRT 
est_teste_frt <- tau/sqrt(n*s2/(n1*n0))
pvalor_frt_clt <- 2 - 2*pnorm(abs(est_teste_frt)) # Bilateral

# FRT Monte Carlo
mc <- 10^4
est_teste_frt_mc <- rep(0, mc)
for (i in 1:mc) {
  zpermut <- sample(z)
  tau_permut <- mean(y[zpermut == 1]) - mean(y[zpermut == 0])
  est_teste_frt_mc[i] <- tau_permut/sqrt(n*s2/(n1*n0))
}
pvalor_frt_mc <- mean(abs(est_teste_frt_mc) > abs(est_teste_frt))

# SRE
z_permuta_SRE <- function(z, x) {
  x_levels <- unique(x)
  K <- length(x_levels)
  z_perm <- z
  for (k in 1:K) {
    x_k <- x_levels[k]
    z_perm[x == x_k] <- sample(z[x == x_k])
  }
  return(z_perm)
}

ee <- function(y, z, x) {
  x_levels <- unique(x)
  K <- length(x_levels)
  pi_k <- rep(0, K)
  tau_k <- rep(0, K)
  n <- length(z)
  for (k in 1:K) {
    x_k <- x_levels[k]
    z_k <- z[x == x_k]
    y_k <- y[x == x_k]
    pi_k[k] <- length(z_k)/n
    tau_k[k] <- mean(y_k[z_k == 1]) - mean(y_k[z_k == 0])
  }
  tau_S <- sum(pi_k * tau_k)
  return(tau_S)
}

# race
mc <- 10^4
estatisticas_x1 <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x1)
  estatisticas_x1[i, 1] <- ee(y, z_perm, x1)
}

# Teste Bilateral 
ee_obs_x1 <- ee(y, z, x1)
pvalor_frt_sre_x1 <- 2 - 2*pnorm(abs(ee_obs_x1))
pvalor_frt_sre_x1_mc <- mean(abs(estatisticas_x1[, 1]) > abs(ee_obs_x1))


# married
x2 <- lalonde$married

estatisticas_x2 <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x2)
  estatisticas_x2[i, 1] <- ee(y, z_perm, x2)
}

# Teste Bilateral 
ee_obs_x2 <- ee(y, z, x2)
pvalor_frt_sre_x2 <- 2 - 2*pnorm(abs(ee_obs_x2))
pvalor_frt_sre_x2_mc <- mean(abs(estatisticas_x2[, 1]) > abs(ee_obs_x2))

# nodegr
x3 <- lalonde$nodegr

estatisticas_x3 <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x3)
  estatisticas_x3[i, 1] <- ee(y, z_perm, x3)
}

# Teste Bilateral 
ee_obs_x3 <- ee(y, z, x3)
pvalor_frt_sre_x3 <- 2 - 2*pnorm(abs(ee_obs_x3))
pvalor_frt_sre_x3_mc <- mean(abs(estatisticas_x3[, 1]) > abs(ee_obs_x3))
```

## Neyman

``` r
# CRE
tau_hat <- mean(y[z == 1]) - mean(y[z == 0])
V_hat <- var(y[z == 1]) / n1 + var(y[z == 0]) / n0
se_hat <- sqrt(V_hat)

est_neyman_cre <- tau_hat/se_hat
ic_neyman_cre <- c(est_neyman_cre - 1.96 * se_hat, est_neyman_cre + 1.96 * se_hat)

# SRE
neyman_sre <- function(y, z, x) {
  x_levels <- unique(x)
  K <- length(x_levels)
  pi_k <- rep(0, K)
  tau_k <- rep(0, K)
  v_k <- rep(0, K)
  n <- length(z)
  for (k in 1:K) {
    x_k <- x_levels[k]
    z_k <- z[x == x_k]
    y_k <- y[x == x_k]
    pi_k[k] <- length(z_k)/n
    tau_k[k] <- mean(y_k[z_k == 1]) - mean(y_k[z_k == 0])
    v_k[k] <- var(y_k[z_k == 1]) / sum(z_k) + var(y_k[z_k == 0]) / sum(1 - z_k)
  }
  tau_s <- sum(pi_k * tau_k)
  v_s <- sum(pi_k^2 * v_k) 
   
  return(c(tau_s, v_s))
}

# race
est_x1 <- neyman_sre(y, z, x1)

se_hat_x1 <- sqrt(est_x1[2])
est_neyman_x1 <- est_x1[1]/se_hat_x1

# IC 95% (convervador)
ic_neyman_sre_x1 <- c(est_neyman_x1 - 1.96 * se_hat_x1, est_neyman_x1 + 1.96 * se_hat_x1)

# married
est_x2 <- neyman_sre(y, z, x2)

se_hat_x2 <- sqrt(est_x2[2])
est_neyman_x2 <- est_x2[1]/se_hat_x2

# IC 95% (convervador)
ic_neyman_sre_x2 <- c(est_neyman_x2 - 1.96 * se_hat_x2, est_neyman_x2 + 1.96 * se_hat_x2)

# nodegr
est_x3 <- neyman_sre(y, z, x3)

se_hat_x3 <- sqrt(est_x3[2])
est_neyman_x3 <- est_x3[1]/se_hat_x3

# IC 95% (convervador)
ic_neyman_sre_x3 <- c(est_neyman_x3 - 1.96 * se_hat_x3, est_neyman_x3 + 1.96 * se_hat_x3)
```

## Resultados

### *race*

``` r
# p-valores FRT
list("FRT CRE" = pvalor_frt_clt,
     "FRT CRE (MC)" = pvalor_frt_mc,
     "FRT SRE (MC)" = pvalor_frt_sre_x1_mc)
```

    ## $`FRT CRE`
    ## [1] 0.00490649
    ## 
    ## $`FRT CRE (MC)`
    ## [1] 0.0048
    ## 
    ## $`FRT SRE (MC)`
    ## [1] 0.0042

``` r
# ICs 95% conservadores Neyman
list("Neyman CRE" = ic_neyman_cre,
     "Neyman SRE" = ic_neyman_sre_x1)
```

    ## $`Neyman CRE`
    ## [1] -1312.479  1317.828
    ## 
    ## $`Neyman SRE`
    ## [1] -1322.155  1327.466

### *marital status*

``` r
# p-valores FRT
list("FRT CRE" = pvalor_frt_clt,
     "FRT CRE (MC)" = pvalor_frt_mc,
     "FRT SRE (MC)" = pvalor_frt_sre_x2_mc)
```

    ## $`FRT CRE`
    ## [1] 0.00490649
    ## 
    ## $`FRT CRE (MC)`
    ## [1] 0.0048
    ## 
    ## $`FRT SRE (MC)`
    ## [1] 0.0057

``` r
# ICs 95% conservadores Neyman
list("Neyman CRE" = ic_neyman_cre,
     "Neyman SRE" = ic_neyman_sre_x2)
```

    ## $`Neyman CRE`
    ## [1] -1312.479  1317.828
    ## 
    ## $`Neyman SRE`
    ## [1] -1311.143  1316.416

### *high school diploma*

``` r
# p-valores FRT
list("FRT CRE" = pvalor_frt_clt,
     "FRT CRE (MC)" = pvalor_frt_mc,
     "FRT SRE (MC)" = pvalor_frt_sre_x3_mc)
```

    ## $`FRT CRE`
    ## [1] 0.00490649
    ## 
    ## $`FRT CRE (MC)`
    ## [1] 0.0048
    ## 
    ## $`FRT SRE (MC)`
    ## [1] 0.0115

``` r
# ICs 95% conservadores Neyman
list("Neyman CRE" = ic_neyman_cre,
     "Neyman SRE" = ic_neyman_sre_x3)
```

    ## $`Neyman CRE`
    ## [1] -1312.479  1317.828
    ## 
    ## $`Neyman SRE`
    ## [1] -1305.000  1309.792

# Questão 6

# Importando dados

``` r
# Bases de dados
dados <- read.table("https://dataverse.harvard.edu/api/access/datafile/7440282")

# Respostas
y <- log(dados$duration)
# Tratamentos
z <- dados$treatment
# Extratos
k <- dados$quarter
k_levels <- unique(k) # Níveis dos extratos

# Matriz de covariáveis
x <- dados[, -c(1, 2, 10)]
x <- scale(x)

# Tamanho total
n <- nrow(dados)
```

## Estratificação

``` r
# FRT
ee <- function(y, z, x) {
  
  x_levels <- unique(x)
  K <- length(x_levels)
  pi_k <- rep(0, K)
  tau_k <- rep(0, K)
  var_k <- rep(0, K)
  n <- length(z)
  for (k in 1:K) {
    x_k <- x_levels[k]
    z_k <- z[x == x_k]
    y_k <- y[x == x_k]
    n_k <- length(z_k)
    n_k1 <- sum(z_k)
    n_k0 <- n_k - n_k1
    pi_k[k] <- n_k/n
    tau_k[k] <- mean(y_k[z_k == 1]) - mean(y_k[z_k == 0])
    s2_k <- var(y_k)
    var_k[k] <- n_k * s2_k/(n_k1*n_k0)
  }
  
  tau_S <- sum(pi_k * tau_k)
  var_S <- sum(pi_k^2 * var_k)
  
  return(c(tau_S, var_S))
  
}

# IC 95% FRT
ee_obs <- ee(y, z, k)
lim_inf <- ee_obs[1] - 1.96*sqrt(ee_obs[2])
lim_sup <- ee_obs[1] + 1.96*sqrt(ee_obs[2])

ic_frt <- c(lim_inf, lim_sup)

# Neyman
neyman_sre <- function(y, z, x) {
  
  x_levels <- unique(x)
  K <- length(x_levels)
  pi_k <- rep(0, K)
  tau_k <- rep(0, K)
  v_k <- rep(0, K)
  n <- length(z)
  
  for (k in 1:K) {
    x_k <- x_levels[k]
    z_k <- z[x == x_k]
    y_k <- y[x == x_k]
    pi_k[k] <- length(z_k)/n
    
    modelo_k <- lm(y_k ~ z_k)
    
    tau_k[k] <- coef(modelo_k)[2]
    v_k[k] <- car::hccm(modelo_k, type = "hc2")[2, 2]
    
  }
  
  tau_s <- sum(pi_k * tau_k)
  v_s <- sum(pi_k^2 * v_k) 
  
  return(c(tau_s, v_s))
  
}

# IC 95% Neyman (convervador)
est_neyman <- neyman_sre(y, z, k)
lim_inf <- est_neyman[1] - 1.96*sqrt(est_neyman[2])
lim_sup <- est_neyman[1] + 1.96*sqrt(est_neyman[2])
ic_neyman <- c(lim_inf, lim_sup)
```

# Ajuste de Regressão

``` r
# Qtd. de extratos
K <- length(k_levels)

# Proporções do tamanho dos extrato
pi_k <- rep(0, K)

# Efeitos causais médios dos extratos
tau_k <- rep(0, K)

# Variância dentro dos extratos
var_k <- rep(0, K)
  
# Ajuste de regressação por extrato
for (estrato in sort(k_levels)) {
  x_k <- x[k == estrato, ]
  z_k <- z[k == estrato]
  y_k <- y[k == estrato]
  pi_k[estrato] <- length(z_k)/n
    
  modelo_k <- lm(y_k ~ z_k*x_k)
    
  tau_k[estrato] <- coef(modelo_k)[2]
  var_k[estrato] <- car::hccm(modelo_k, type = "hc2")[2, 2]
    
}
  
# Efeito causal médio SRE
tau_S <- sum(pi_k * tau_k)
se_S <- sqrt(sum(pi_k^2 * var_k))

# IC 95% ReM
lim_inf <- tau_S - 1.96 * se_S
lim_sup <- tau_S + 1.96 * se_S
ic_rem <- c(lim_inf, lim_sup)
```

## Resultados

``` r
# Intervalos de confiança 95% (conservador)
list("FRT SRE" = ic_frt,
     "Neyman SRE" = ic_neyman,
     "Lin SRE", ic_rem)
```

    ## $`FRT SRE`
    ## [1] -0.15031236 -0.02950056
    ## 
    ## $`Neyman SRE`
    ## [1] -0.15027005 -0.02954287
    ## 
    ## [[3]]
    ## [1] "Lin SRE"
    ## 
    ## [[4]]
    ## [1] -0.14922947 -0.03129349
