# Gaussian Process implementation details

## Overview

We make use of Gaussian Processes in several places in `EpiNow2`. For
example, the default model for
[`estimate_infections()`](https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.md)
uses a Gaussian Process to model the 1st order difference on the log
scale of the reproduction number. This vignette describes the
implementation details of the approximate Gaussian Process used in
`EpiNow2`.

## Definition

The one-dimensional Gaussian Processes ($`\mathrm{GP}_t`$) we use can be
written as

``` math
\begin{equation}
\mathrm{GP}(\mu(t), k(t, t'))
\end{equation}
```

where $`\mu(t)`$ and $`k(t,t')`$ are the mean and covariance functions,
respectively. In our case as set out above, we have

``` math
\begin{align}
\mu(t) &\equiv 0 \\
k(t,t') &= k(|t - t'|) = k(\Delta t)
\end{align}
```

with the following choices available for the kernel $`k`$

### Matérn 3/2 covariance kernel (the default)

``` math
\begin{equation}
k(\Delta t) = \alpha^2 \left( 1 + \frac{\sqrt{3} \Delta t}{\rho} \right) \exp \left( - \frac{\sqrt{3} \Delta t}{\rho}\right)
\end{equation}
```

with $`\rho>0`$ and $`\alpha > 0`$ the length scale and magnitude,
respectively, of the kernel. Note that here and later we use a slightly
different definition of $`\alpha`$ compared to Riutort-Mayol et
al.^(\[1\]), where this is defined as our $`\alpha^2`$.

### Squared exponential kernel

``` math
\begin{equation}
k(\Delta t) = \alpha^2 \exp \left( - \frac{1}{2} \frac{(\Delta t^2)}{l^2} \right)
\end{equation}
```

### Ornstein-Uhlenbeck (Matérn 1/2) kernel

``` math
\begin{equation}
k(\Delta t) = \alpha^2 \exp{\left( - \frac{\Delta t}{2 \rho^2}  \right)}
\end{equation}
```

### Matérn 5/2 covariance kernel

``` math
\begin{equation}
k(\Delta t) = \alpha \left( 1 + \frac{\sqrt{5} \Delta t}{\rho} + \frac{5}{3} \left(\frac{\Delta t}{l} \right)^2 \right) \exp \left( - \frac{\sqrt{5} \Delta t}{\rho}\right)
\end{equation}
```

## Hilbert space approximation

In order to make our models computationally tractable, we approximate
the Gaussian Process using a Hilbert space approximation to the Gaussian
Process^(\[1\]), centered around mean zero.

``` math
\begin{equation}
\mathrm{GP}(0, k(\Delta t)) \approx \sum_{j=1}^m \left(S_k(\sqrt{\lambda_j}) \right)^\frac{1}{2} \phi_j(t) \beta_j
\end{equation}
```

with $`m`$ the number of basis functions to use in the approximation,
which we calculate from the number of time points $`t_\mathrm{GP}`$ to
which the Gaussian Process is being applied (rounded up to give an
integer value), as is recommended^(\[1\]).

``` math
\begin{equation}
m = b t_\mathrm{GP}
\end{equation}
```

and values of $`\lambda_j`$ given by

``` math
\begin{equation}
\lambda_j = \left( \frac{j \pi}{2 L} \right)^2
\end{equation}
```

where $`L`$ is a positive number termed boundary condition, and
$`\beta_{j}`$ are regression weights with standard normal prior

``` math
\begin{equation}
\beta_j \sim \mathrm{Normal}(0, 1)
\end{equation}
```

The function $`S_k(x)`$ is the spectral density relating to a particular
covariance function $`k`$. In the case of the Matérn kernel of order
$`\nu`$ this is given by

``` math
\begin{equation}
S_k(x) = \alpha^2 \frac{2 \sqrt{\pi} \Gamma(\nu + 1/2) (2\nu)^\nu}{\Gamma(\nu) \rho^{2 \nu}} \left( \frac{2 \nu}{\rho^2} + \omega^2 \right)^{-\left( \nu + \frac{1}{2}\right )}
\end{equation}
```

For $`\nu = 3 / 2`$ (the default in `EpiNow2`) this simplifies to

``` math
\begin{equation}
    S_k(\omega) =
    \alpha^2 \frac{4 \left(\sqrt{3} / \rho \right)^3}{\left(\left(\sqrt{3} / \rho\right)^2 + \omega^2\right)^2} =
    \left(\frac{2 \alpha \left(\sqrt{3} / \rho \right)^{\frac{3}{2}}}{\left(\sqrt{3} / \rho\right)^2 + \omega^2} \right)^2
\end{equation}
```

For $`\nu = 1 / 2`$ it is

``` math
\begin{equation}
S_k(\omega) = \alpha^2 \frac{2}{\rho \left(1 / \rho^2 + \omega^2\right)}
\end{equation}
```

and for $`\nu = 5 / 2`$ it is

``` math
\begin{equation}
S_k(\omega) =
    \alpha^2 \frac{16 \left(\sqrt{5} / \rho \right)^5}{3 \left(\left(\sqrt{5} / \rho\right)^2 + \omega^2\right)^3}
\end{equation}
```

In the case of a squared exponential the spectral kernel density is
given by

``` math
\begin{equation}
S_k(\omega) = \alpha^2 \sqrt{2\pi} \rho \exp \left( -\frac{1}{2} \rho^2 \omega^2 \right)
\end{equation}
```

The functions $`\phi_{j}(x)`$ are the eigenfunctions of the Laplace
operator,

``` math
\begin{equation}
\phi_j(t) = \frac{1}{\sqrt{L}} \sin\left(\sqrt{\lambda_j} (t^* + L)\right)
\end{equation}
```

with time rescaled linearly to be between -1 and 1,

``` math
\begin{equation}
t^* = \frac{t - \frac{1}{2}t_\mathrm{GP}}{\frac{1}{2}t_\mathrm{GP}}
\end{equation}
```

Relevant default priors are

``` math
\begin{align}
\alpha &\sim \mathrm{HalfNormal}(0, 0.01) \\
\rho   &\sim \mathrm{LogNormal} (\mu_\rho, \sigma_\rho)\\
\end{align}
```

with $`\rho`$ additionally constrained with an upper bound of $`60`$ and
$`\mu_{\rho}`$ and $`\sigma_\rho`$ calculated using a mean of 21 and
standard deviation of 7.

Furthermore, by default we set.

``` math
\begin{align}
b &= 0.2 \\
L &= 1.5
\end{align}
```

These values as well as the prior distributions of relevant parameters
can all be changed by the user.

## References

1\.

Riutort-Mayol, G., Bürkner, P.-C., Andersen, M. R., Solin, A., &
Vehtari, A. (2020). *Practical hilbert space approximate bayesian
gaussian processes for probabilistic programming*.
<https://arxiv.org/abs/2004.11408>
