# Asymptotic covariance for the PMM3 regression estimator

Working notes that derive the asymptotic covariance matrix used by
`vcov.PMM3fit()`. The result mirrors the PMM2 case of Zabolotnii et al.
(2018, 2022, 2023) and specialises Kunchenko's polynomial-maximisation
framework to symmetric platykurtic errors (s = 3, with all odd-order
moments equal to zero).

## Setup

Consider the linear model
$$ y_i = \mathbf{x}_i^\top \boldsymbol\beta + \varepsilon_i,
   \qquad i = 1, \dots, n, $$
with i.i.d. errors $\varepsilon_i$ satisfying
$\mathrm{E}[\varepsilon_i] = 0$,
$\mathrm{E}[\varepsilon_i^{2k+1}] = 0$ (symmetry),
$\mathrm{Var}(\varepsilon_i) = \sigma^2 = m_2$,
$m_4 = \mathrm{E}[\varepsilon_i^4]$,
$m_6 = \mathrm{E}[\varepsilon_i^6]$.

The standardised cumulant coefficients are
$$ \gamma_4 = \frac{m_4}{m_2^2} - 3, \qquad
   \gamma_6 = \frac{m_6}{m_2^3} - 15\,\frac{m_4}{m_2^2} + 30. $$

PMM3 uses the third-order polynomial score
$$ S_i(\boldsymbol\beta) =
   a_1 (y_i - \mathbf{x}_i^\top\boldsymbol\beta)
 + a_3 (y_i - \mathbf{x}_i^\top\boldsymbol\beta)^3, $$
with the optimal coefficients $a_1, a_3$ chosen to minimise the
asymptotic variance of the resulting estimator within the class of
symmetric polynomial estimators of degree $\le 3$. Asymmetric terms
($k = 2$) are excluded by the symmetry assumption.

## Derivation outline

Following Kunchenko's classical argument (see Zabolotnii et al. 2018,
Theorem 1, and 2022, §3) the optimal coefficients are obtained from the
moment system
$$ \mathrm{E}\!\left[S_i(\boldsymbol\beta)\,\dot S_i(\boldsymbol\beta)\right]
   = - \nabla_{\boldsymbol\beta}\, \mathrm{E}[S_i(\boldsymbol\beta)], $$
which for symmetric errors reduces to a 2x2 linear system in
$(a_1, a_3)$ whose matrix has entries
$$ \begin{pmatrix} m_2 & 3 m_4 \\ 3 m_4 & 15 m_6 \end{pmatrix}. $$
Solving this system and substituting into the asymptotic variance
formula
$$ \mathrm{Avar}(\hat{\boldsymbol\beta}_{\mathrm{PMM3}}) \;=\;
   \frac{m_2}{\mathrm{Var}[S_i] / m_2^2}\,(X^\top X)^{-1} $$
yields, after simplification,
$$ \boxed{\;
   \mathrm{Avar}(\hat{\boldsymbol\beta}_{\mathrm{PMM3}})
   \;=\; g_3 \, \sigma^2 \, (X^\top X)^{-1},
   \qquad
   g_3 \;=\; 1 \;-\; \frac{\gamma_4^{\,2}}{6 + 9\gamma_4 + \gamma_6}
   \;} $$
provided that the denominator $6 + 9\gamma_4 + \gamma_6 > 0$, which is
the standard PMM-feasibility condition for symmetric distributions.

The structure is identical to PMM2 -- only the scalar efficiency
factor changes:

| Estimator | Efficiency factor                                            |
|-----------|--------------------------------------------------------------|
| OLS       | $g = 1$                                                      |
| PMM2      | $g_2 = 1 - c_3^{2} / (2 + c_4)$                              |
| PMM3      | $g_3 = 1 - \gamma_4^{2} / (6 + 9\gamma_4 + \gamma_6)$        |

Under Gaussian errors $\gamma_4 = 0$, hence $g_3 = 1$ and PMM3 reduces
to OLS asymptotically.

## Practical implications

* For symmetric platykurtic distributions ($\gamma_4 < 0$), PMM3 gives
  $g_3 < 1$, i.e. strictly smaller asymptotic variance than OLS,
  provided the feasibility condition holds.
* For symmetric mesokurtic ($\gamma_4 \approx 0$) errors there is no
  efficiency gain over OLS.
* PMM3 is **not** robust to asymmetry: when $\gamma_3 \neq 0$, the
  estimator remains consistent but the variance formula above is
  invalid. The package therefore emits a warning whenever
  $|\hat\gamma_3| > 0.3$ and recommends PMM2.

## Implementation

The package implements the formula above in
`pmm3_variance_matrices(X, m2, m4, m6)`, returning both the OLS
covariance $m_2 (X^\top X)^{-1}$ and the PMM3 covariance
$g_3 m_2 (X^\top X)^{-1}$. The S4 method `vcov(PMM3fit)` calls this
helper using the stored model matrix and the PMM3 moments computed on
the initial OLS residuals.

## References

* Kunchenko, Yu. P. *Polynomial Parameter Estimations of Close to
  Gaussian Random Variables*. Kyiv, 2002.
* Zabolotnii, S. W., Warsza, Z. L., & Tkachenko, O. (2018). On
  Estimation of Linear Regression Parameters... *Advances in
  Intelligent Systems and Computing*.
  doi:10.1007/978-3-319-77179-3_75
* Zabolotnii, S. W., Tkachenko, O., & Warsza, Z. L. (2022).
  doi:10.1007/978-3-031-03502-9_37
* Zabolotnii, S. W., Tkachenko, O., & Warsza, Z. L. (2023).
  doi:10.1007/978-3-031-25844-2_21
