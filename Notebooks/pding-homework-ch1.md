Ding, P. (2023) 1.5 Homework Problems
================

Homework problems from [Ding, P. (2023)](#1), Chapter 1: Correlation,
Assossiation, and the Yule-Simpson Paradox.

### *1.1 Independence in two-by-two tables*

Prove (1) and (2) in Proposition 1.1.

1)  The following statements are equivalent: $Z \perp \!\!\! \perp Y$,
    $RD = 0$, $RR = 1$ and $OR = 1$.

$$Z \perp \!\!\! \perp Y \Rightarrow \textrm{pr}(Y|Z) = \textrm{pr}(Y)$$
$$RD = pr(Y = 1 | Z = 1) - pr(Y = 1 |Z = 0) \overset{Z \perp \!\!\! \perp Y}{=} pr(Y = 1) - pr(Y = 1) = 0$$

$$RR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} \overset{Z \perp \!\!\! \perp Y}{=} = \frac{pr(Y = 1)}{pr(Y = 1) = 1}$$

$$OR = \frac{\frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)}}{\frac{pr(Y = 0 | Z = 1)}{pr(Y = 0 |Z = 0)}} \overset{Z \perp \!\!\! \perp Y}{=} \frac{\frac{pr(Y = 1)}{pr(Y = 1)}}{\frac{pr(Y = 0)}{pr(Y = 0)}} = 1$$

2)  If $p_{zy}$’s are all positive, then $RD > 0$ is equivalent to
    $RR > 1$ and is also equivalent to $OR > 1$.

$$p_{zy} > 0,\; RD > 0 \Rightarrow pr(Y = 1|Z = 1) > pr(Y = 1|Z = 0)$$
$$RR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} > 1$$
$$pr(Y = 1|Z = 1) > pr(Y = 1|Z = 0) \Rightarrow \frac{p_{11}}{pr(Z=1)} > \frac{p_{10}}{pr(Z=0)} \Rightarrow \frac{-p_{10}}{pr(Z=1)} > \frac{-p_{00}}{pr(Z=0)} \Rightarrow pr(Y = 0|Z = 1) < pr(Y = 0|Z = 0)$$
$$OR = \frac{pr(Y = 1 | Z = 1)}{pr(Y = 1 |Z = 0)} \cdot \frac{pr(Y = 0 | Z = 0)}{pr(Y = 0 |Z = 1)} >1$$
\### References

\[1\] <a id="1">Ding, P. (2023)</a>. **A First Course in Causal
Inference**, Chapman and Hall/CRC. (<https://arxiv.org/abs/2305.18793>)
