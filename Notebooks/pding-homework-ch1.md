Ding, P. (2023) 1.5 Homework Problems
================

Homework problems from [Ding, P. (2023)](#1), Chapter 1: Correlation,
Assossiation, and the Yule-Simpson Paradox.

### *1.1 Independence in two-by-two tables*

Prove (1) and (2) in Proposition 1.1.

**Solution**

1)  The following statements are equivalent: $Z \perp Y$, $RD = 0$,
    $RR = 1$ and $OR = 1$.

$$Z \perp \!\!\! \perp Y \Rightarrow \textrm{pr}(Y|Z) = \textrm{pr}(Y)$$
$$RD = pr(Y = 1 | Z = 1) - pr(Y = 1 |Z = 0) = pr(Y = 1) - pr(Y = 1) = 0$$

$$RR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} = \frac{pr(Y = 1)}{pr(Y = 1) = 1}$$

$$OR = \frac{\frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)}}{\frac{pr(Y = 0 | Z = 1)}{pr(Y = 0 |Z = 0)}} = \frac{\frac{pr(Y = 1)}{pr(Y = 1)}}{\frac{pr(Y = 0)}{pr(Y = 0)}} = 1$$

2)  If $p_{zy}$’s are all positive, then $RD > 0$ is equivalent to
    $RR > 1$ and is also equivalent to $OR > 1$

$$p_{zy} > 0,\; RD > 0 \Rightarrow pr(Y = 1|Z = 1) > pr(Y = 1|Z = 0)$$
$$RR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} > 1$$
$$pr(Y = 1|Z = 1) > pr(Y = 1|Z = 0) \Rightarrow \frac{p_{11}}{pr(Z=1)} > \frac{p_{10}}{pr(Z=0)} \Rightarrow \frac{-p_{10}}{pr(Z=1)} > \frac{-p_{00}}{pr(Z=0)} \Rightarrow pr(Y = 0|Z = 1) < pr(Y = 0|Z = 0)$$
$$OR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} \cdot \frac{pr(Y = 0 | Z = 0)}{pr(Y = 0 |Z = 1)} >1$$

### *1.2 More examples of the Yule–Simpson Paradox*

Give a numeric example of a two-by-two-by-two table in which the Yule–
Simpson Paradox arises. Find a real-life example in which the
Yule–Simpson Paradox arises

**Solution**

For a numeric example, assume that the two-by-two table based on the
aggregated data has counts

| whole population | $Y = 1$ | $Y = 0$ |
|:----------------:|:-------:|:-------:|
|     $Z = 1$      |   30    |   20    |
|     $Z = 0$      |   40    |   10    |

$$\widehat{RD} = \frac{30}{50} - \frac{40}{50} = 60\\% - 80\\% = -20\\% < 0$$

The two two-by-two table based on the subgroup with $X = 1$ has counts

| subpopulation $X = 1$ | $Y = 1$ | $Y = 0$ |
|:---------------------:|:-------:|:-------:|
|        $Z = 1$        |   18    |   11    |
|        $Z = 0$        |   12    |    9    |

$$\widehat{RD}_{X=1} = \frac{18}{29} - \frac{12}{21} \approx 62\\% - 57\\% = 5\\% > 0$$
For the subgroup with $X = 0$, the two two-by-two table has counts

| subpopulation $X = 0$ | $Y = 1$ | $Y = 0$ |
|:---------------------:|:-------:|:-------:|
|        $Z = 1$        |   20    |   12    |
|        $Z = 0$        |   10    |    8    |

$$\widehat{RD}_{X=0} = \frac{20}{32} - \frac{10}{18} \approx 62\\% - 56\\% = 6\\% > 0$$
A real-life example of Simpson’s paradox is discussed in [Wang, Z.,
Rousseau, R. (2021)](#2), where the infection fatality rate of COVID-19
in men versus women is debated. In summary, the analysis of aggregated
data pointed to a higher infection fatality rate for women. However,
when evaluated by age group, the fatality rate was higher for men across
all age groups.

### *1.3 Correlation and partial correlation*

Consider a three-dimensional Normal random vector:

$$\begin{pmatrix}X\\ Y\\ Z\end{pmatrix} \sim N\begin{pmatrix}\begin{pmatrix}0\\ 0\\ 0\end{pmatrix},\begin{pmatrix}1 &\rho_{XY} & \rho_{XZ}\\ \rho_{XY} &  1& \rho_{YZ}\\ \rho_{XZ} & \rho_{YZ} & 1\end{pmatrix}\end{pmatrix}$$
The correlation coefficient between Y and Z is ρY Z. There are many
equivalent definitions of the partial correlation coefficient. For a
multivariate Normal vector, let ρY Z\|X denote the partial correlation
coefficient between Y and Z given X, which is defined as their
correlation coefficient in the conditional distribution (Y, Z) \| X.
Show that

$$\rho_{YZ|X} = \frac{\rho_{YZ} - \rho_{YX}\rho_{ZX}}{\sqrt{(1-\rho^{2}_{YX})(1-\rho^{2}_{ZX})}}$$
Give a numerical example with $\rho_{YZ} > 0$ and $\rho_{YZ|X} < 0$.

Remark: This is the Yule–Simpson Paradox for a Normal random vector. You
can use the results in Chapter A.1.2 to prove the formula for the
partial correlation coefficient.

**Solution**

### References

\[1\] <a id="2">Ding, P. (2023)</a>. **A First Course in Causal
Inference**, Chapman and Hall/CRC. (<https://arxiv.org/abs/2305.18793>)

\[2\] <a id="1">Wang, Z., Rousseau, R. (2021)</a>. **COVID-19, the
Yule-Simpson paradox and research evaluation**. Scientometrics 126,
3501–3511 (2021). (<https://doi.org/10.1007/s11192-020-03830-w>)
