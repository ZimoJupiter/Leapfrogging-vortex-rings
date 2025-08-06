# Simulations for Leapfrogging vortex rings using Bito-Savart law and Lamb-Ossen vortex core model.

## Biot-Savart Law

The Biot-Savart Law in fluid mechanics emerges naturally from the solution of the Poisson equation for vorticity-induced velocity fields. Applying the vector identity $\nabla \times (\nabla \times \mathbf{V}) = \nabla(\nabla \cdot \mathbf{V}) - \nabla^2 \mathbf{V}$ and the incompressibility condition yields the vector Poisson equation for velocity:

$$
\nabla^2 \mathbf{V} = -\nabla \times \boldsymbol{\omega}
$$

The fundamental solution to the Poisson equation in three-dimensional free space is given by the Green's function:

$$
G(\mathbf{r}) = \frac{1}{4\pi|\mathbf{r}|}, \quad \text{satisfying} \quad \nabla^2 G(\mathbf{r}) = -\delta(\mathbf{r})
$$

The velocity field solution can then be expressed as:

$$
\mathbf{V}(\mathbf{r}) = \int_{\mathbb{R}^3} \nabla \times \boldsymbol{\omega}(\mathbf{r}') G(\mathbf{r}-\mathbf{r}') \, \text{d}V'
$$

where $\mathbf{r}'$ and $\mathbf{r}$ represent the position vectors of the source point and the observation point, respectively. This leads to the classical form:

$$
\mathbf{V}(\mathbf{r}) = \frac{1}{4\pi} \int_{\mathbb{R}^3} \frac{\boldsymbol{\omega}(\mathbf{r}') \times (\mathbf{r} - \mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|^3} \, \text{d}V'
$$

For concentrated vortex lines (where vorticity $\boldsymbol{\omega}$ is confined to a curve $C$ with circulation $\Gamma$), this reduces to:

$$
\mathbf{V}(\mathbf{r}) = \frac{\Gamma}{4\pi} \oint_C \frac{\text{d}\mathbf{l}' \times (\mathbf{r} - \mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|^3}
$$

![image](https://github.com/ZimoJupiter/Leapfrogging-vortex-rings/blob/main/Plots/BS%20law.png)

## Vortex Core Correction

To avoid singularities near the vortex line, vortex core model corrections are adopted:

$$
\mathbf{V}(\mathbf{r}) = \frac{\Gamma}{4\pi} \oint_C \frac{h^2}{\left(r_c^{2n}+h^{2n}\right)^{1/n}} \cdot \frac{\text{d}\mathbf{l}' \times (\mathbf{r} - \mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|^3}
$$

where:
- $h$ = distance between point of interest and filaments
- $n$ = vortex core type parameter (Lamb-Oseen uses $n=2$)
- $r_c$ = vortex core radius, evolving as:

$$
r_c(\zeta) = \sqrt{r_0^2\frac{|\mathbf{l}|}{\Delta |\mathbf{l}| + |\mathbf{l}|} + \frac{4\alpha_L\delta\nu\zeta}{\Omega}}
$$

Empirical parameters:
- $\delta$ ≈ 10-1000 (viscosity coefficient)
- $r_0$ ≈ 5-10% of local chord length (initial core radius)

![image](https://github.com/ZimoJupiter/Leapfrogging-vortex-rings/blob/main/Plots/Vortex%20core.png)

## Leapfrogging vortex rings

The vortex ring positions are updated using an explicit algorithm. The animation is as follows:

![image](https://github.com/ZimoJupiter/Leapfrogging-vortex-rings/blob/main/Plots/Animation.gif)
