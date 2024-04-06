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

``` r
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
n_controle <- sum(df_multicenter$n0)
n_1mg <- sum(df_multicenter$n1)
n_5mg <- sum(df_multicenter$n5)

total_1mg <- n_controle + n_1mg

df_multicenter_1mg <- df_multicenter %>% 
  select(colnames(df_multicenter)[!str_detect(colnames(df_multicenter), "5$")]) %>% 
  mutate(pi_k = (n0+n1)/total_1mg,
         tau_k = mean1 - mean0,
         var_k = sd1^2 / n1 + sd0^2 / n0)

tau_s_1mg <- sum(df_multicenter_1mg$pi_k * df_multicenter_1mg$tau_k)
v_s_1mg <- sum(df_multicenter_1mg$pi_k^2 * df_multicenter_1mg$var_k)

lim_inf_1mg <- tau_s_1mg - 1.96*sqrt(v_s_1mg)
lim_sup_1mg <- tau_s_1mg + 1.96*sqrt(v_s_1mg)
c(lim_inf_1mg, lim_sup_1mg)
```

    ## [1] -1.468874  0.154194

``` r
est_s_1mg <- tau_s_1mg/sqrt(v_s_1mg)

2 - 2 * pnorm(abs(est_s_1mg))
```

    ## [1] 0.1123782

``` r
total_5mg <- n_controle + n_5mg

df_multicenter_5mg <- df_multicenter %>% 
  select(colnames(df_multicenter)[!str_detect(colnames(df_multicenter), "1$")]) %>% 
  mutate(pi_k = (n0+n5)/total_5mg,
         tau_k = mean5 - mean0,
         var_k = sd5^2 / n5 + sd0^2 / n0)

tau_s_5mg <- sum(df_multicenter_5mg$pi_k * df_multicenter_5mg$tau_k)
v_s_5mg <- sum(df_multicenter_5mg$pi_k^2 * df_multicenter_5mg$var_k)

lim_inf_5mg <- tau_s_5mg - 1.96*sqrt(v_s_5mg)
lim_sup_5mg <- tau_s_5mg + 1.96*sqrt(v_s_5mg)
c(lim_inf_5mg, lim_sup_5mg)
```

    ## [1] -2.4259519 -0.8555447

``` r
est_s_5mg <- tau_s_5mg/sqrt(v_s_5mg)

2 - 2 * pnorm(abs(est_s_5mg))
```

    ## [1] 4.211066e-05

# Questão 3

``` r
data("lalonde")

lalonde <- lalonde %>% 
  mutate(race = ifelse(black == 1, "black",
                       ifelse(hisp == 1, "hisp", "other")))

z <- lalonde$treat
x1 <- lalonde$race
y <- lalonde$re78

n <- nrow(lalonde)
n1 <- sum(z)
n0 <- n - n1
tau <- mean(y[z == 1]) - mean(y[z == 0])
s2 = var(y)
```

``` r
# FRT 
est_teste_frt <- tau/sqrt(n*s2/(n1*n0))
pvalor_frt_clt <- 2 - 2*pnorm(est_teste_frt) # Bilateral

# FRT Monte Carlo
mc <- 10^4
est_teste_frt_mc <- rep(0, mc)
for (i in 1:mc) {
  zpermut <- sample(z)
  tau_permut <- mean(y[zpermut == 1]) - mean(y[zpermut == 0])
  est_teste_frt_mc[i] <- tau_permut/sqrt(n*s2/(n1*n0))
}
pvalor_frt_mc <- mean(abs(est_teste_frt_mc) > abs(est_teste_frt))


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

mc <- 10^4
estatisticas <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x1)
  estatisticas[i, 1] <- ee(y, z_perm, x1)
}

ee_obs <- ee(y, z, x1)

mean(estatisticas[, 1] < ee_obs)
```

    ## [1] 0.9976

``` r
x2 <- lalonde$married

estatisticas_x2 <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x2)
  estatisticas_x2[i, 1] <- ee(y, z_perm, x2)
}

ee_obs_x2 <- ee(y, z, x2)

mean(estatisticas_x2[, 1] < ee_obs_x2)
```

    ## [1] 0.9968

``` r
x3 <- lalonde$nodegr

estatisticas_x3 <- matrix(NA, ncol = 2, nrow = mc)
for (i in 1:mc) {
  z_perm <- z_permuta_SRE(z, x3)
  estatisticas_x3[i, 1] <- ee(y, z_perm, x3)
}

ee_obs_x3 <- ee(y, z, x3)

mean(estatisticas_x3[, 1] < ee_obs_x3)
```

    ## [1] 0.9939

``` r
# CRE
tau_hat <- mean(y[z == 1]) - mean(y[z == 0])
V_hat <- var(y[z == 1]) / n1 + var(y[z == 0]) / n0
se_hat <- sqrt(V_hat)

est_neyman <- tau_hat/se_hat

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

est_x1 <- neyman_sre(y, z, x1)

est_x1[1]/sqrt(est_x1[2])
```

    ## [1] 2.65558

``` r
est_x2 <- neyman_sre(y, z, x2)

est_x2[1]/sqrt(est_x2[2])
```

    ## [1] 2.636411

``` r
est_x3 <- neyman_sre(y, z, x3)

est_x3[1]/sqrt(est_x3[2])
```

    ## [1] 2.396085

# Questão 6

``` r
dados <- read.table("https://dataverse.harvard.edu/api/access/datafile/7440282")
```
