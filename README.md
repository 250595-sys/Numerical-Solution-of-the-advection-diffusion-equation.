
# Numerical-Solution-of-the-advection-diffusion-equation.
In this project, I used range of different method to solve advection-diffusion equation. In addition, there is a comparison of all methods, drawing their graphs.

# Numerical Schemes for the 1D Linear Advection Equation

This document provides mathematical formulations, amplification factors, stability analysis, and dissipation/dispersion properties for six classical finite difference schemes solving the 1D linear advection equation:

$$\frac{\partial U}{\partial t} + v \frac{\partial U}{\partial x} = 0$$

where $v > 0$ is the constant wave speed, $\Delta t$ is the time step, $\Delta x$ is the spatial step, and $\sigma = \frac{v \Delta t}{\Delta x}$ is the Courant-Friedrichs-Lewy (CFL) parameter.

---

## 1. FTCS (Forward Time Centered Space)

### Definition

An explicit, first-order in time and second-order in space finite difference scheme. It approximates the temporal derivative using forward differences and the spatial derivative using centered differences.

### Mathematical Formula

$$U_j^{n+1} = U_j^n - \frac{\sigma}{2} \left( U_{j+1}^n - U_{j-1}^n \right)$$

### Amplification Factor ($G$)

Substituting the von Neumann mode $U_j^n = G^n e^{i j \theta}$ (where $\theta = k \Delta x$ is the wave phase angle):

$$G(\theta) = 1 - i \sigma \sin\theta$$

$$\vert{}G(\theta)\vert{}^2 = 1 + \sigma^2 \sin^2\theta > 1 \quad \forall \theta \neq 0$$

### Dissipation and Dispersion

* **Dissipation:** Negative numerical dissipation (unstable amplification). Energy grows exponentially without bound.
* **Dispersion:** Possesses phase error, but it is masked by rapid numerical instability.
* **Stability:** **Unconditionally Unstable**.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Transmissive / Zero-Gradient ($U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$)

---

## 2. Upwind Method (First-Order Upwind)

### Definition

An explicit, first-order space and time scheme that respects the physical direction of wave propagation (information flow) by taking spatial differences from the side from which the fluid/wave arrives.

### Mathematical Formula (for $v > 0$)

$$U_j^{n+1} = U_j^n - \sigma \left( U_j^n - U_{j-1}^n \right)$$

### Amplification Factor ($G$)

$$G(\theta) = (1 - \sigma + \sigma \cos\theta) - i \sigma \sin\theta$$

$$\vert{}G(\theta)\vert{}^2 = 1 - 2\sigma(1 - \sigma)(1 - \cos\theta)$$

### Dissipation and Dispersion

* **Dissipation:** Very high numerical diffusion/viscosity (truncation error generates an artificial diffusion term $\frac{v \Delta x}{2}(1 - \sigma) \frac{\partial^2 U}{\partial x^2}$). Waves undergo heavy amplitude dampening and smoothing over time.
* **Dispersion:** Moderate phase lagging/leading errors depending on $\sigma$. At $\sigma = 1$, dispersion and dissipation vanish completely (exact solution).
* **Stability:** Conditionally stable for $0 \le \sigma \le 1$.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Transmissive / Zero-Gradient ($U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$)

---

## 3. Lax-Friedrichs Method

### Definition

An explicit, first-order scheme created to stabilize the FTCS scheme by replacing the center value $U_j^n$ with the spatial average of its two immediate neighbors $\frac{1}{2}(U_{j+1}^n + U_{j-1}^n)$.

### Mathematical Formula

$$U_j^{n+1} = \frac{1}{2} \left( U_{j+1}^n + U_{j-1}^n \right) - \frac{\sigma}{2} \left( U_{j+1}^n - U_{j-1}^n \right)$$

### Amplification Factor ($G$)

$$G(\theta) = \cos\theta - i \sigma \sin\theta$$

$$\vert{}G(\theta)\vert{}^2 = \cos^2\theta + \sigma^2 \sin^2\theta = 1 - (1 - \sigma^2)\sin^2\theta$$

### Dissipation and Dispersion

* **Dissipation:** Very strong numerical diffusion (adds artificial viscosity proportional to $\frac{\Delta x^2}{2\Delta t}$). Sharp gradients flatten rapidly.
* **Dispersion:** High phase errors for short-wavelength modes.
* **Stability:** Conditionally stable for $\sigma \le 1$.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Transmissive / Zero-Gradient ($U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$)

---

## 4. Lax-Wendroff Method

### Definition

An explicit, second-order accurate in space and time scheme derived from the Taylor series expansion up to second order, using the PDE to replace second time derivatives with spatial derivatives ($\frac{\partial^2 U}{\partial t^2} = v^2 \frac{\partial^2 U}{\partial x^2}$).

### Mathematical Formula

$$U_j^{n+1} = U_j^n - \frac{\sigma}{2}\left(U_{j+1}^n - U_{j-1}^n\right) + \frac{\sigma^2}{2}\left(U_{j+1}^n - 2U_j^n + U_{j-1}^n\right)$$

### Amplification Factor ($G$)

$$G(\theta) = 1 - \sigma^2(1 - \cos\theta) - i \sigma \sin\theta$$

$$\vert{}G(\theta)\vert{}^2 = 1 - \sigma^2(1 - \sigma^2)(1 - \cos\theta)^2$$

### Dissipation and Dispersion

* **Dissipation:** Low dissipation (fourth-order damping term). Preserves peak amplitude much better than Upwind or Lax-Friedrichs.
* **Dispersion:** Predominantly **dispersive error**. Because high-frequency waves travel slower than the analytical phase speed, trailing oscillations (wiggles) appear near sharp fronts.
* **Stability:** Conditionally stable for $\sigma \le 1$.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Transmissive / Zero-Gradient ($U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$)

---

## 5. Leapfrog Method (CTCS - Centered Time Centered Space)

### Definition

An explicit, second-order 3-level scheme that uses a central difference approximation for both time and space derivatives across three temporal states ($n+1$, $n$, and $n-1$).

### Mathematical Formula

$$U_j^{n+1} = U_j^{n-1} - \sigma \left( U_{j+1}^n - U_{j-1}^n \right)$$

### Amplification Factor ($G$)

Solving the quadratic characteristic equation $G^2 + 2i\sigma\sin\theta\,G - 1 = 0$:

$$G(\theta) = -i \sigma \sin\theta \pm \sqrt{1 - \sigma^2 \sin^2\theta}$$

$$\vert{}G(\theta)\vert{} = 1 \quad \text{for } \sigma \le 1$$

### Dissipation and Dispersion

* **Dissipation:** **Zero numerical dissipation** ($\vert{}G\vert{} = 1$). Wave amplitudes do not decay over time.
* **Dispersion:** Significant phase speed errors lead to trailing dispersive ripples. Additionally, because it uses 3 time levels, it supports an unphysical **computational mode** that causes odd/even grid point decoupling.
* **Stability:** Conditionally stable for $\sigma \le 1$. Requires a 2-level scheme (like Upwind) to perform the very first step $U^1$.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Transmissive / Zero-Gradient ($U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$)

---

## 6. BTCS (Backward Time Centered Space - Implicit)

### Definition

An implicit, first-order in time and second-order in space scheme. The spatial difference is evaluated at the new time level $n+1$, resulting in a tridiagonal linear system that must be solved at every time step.

### Mathematical Formula

$$-\frac{\sigma}{2} U_{j-1}^{n+1} + U_j^{n+1} + \frac{\sigma}{2} U_{j+1}^{n+1} = U_j^n$$

In matrix form: $A \mathbf{U}^{n+1} = \mathbf{U}^n$

### Amplification Factor ($G$)

$$G(\theta) = \frac{1}{1 + i \sigma \sin\theta}$$

$$\vert{}G(\theta)\vert{}^2 = \frac{1}{1 + \sigma^2 \sin^2\theta} \le 1 \quad \forall \sigma$$

### Dissipation and Dispersion

* **Dissipation:** Moderate to high numerical dissipation, especially at large Courant numbers ($\sigma > 1$).
* **Dispersion:** Significant phase lag at larger time steps.
* **Stability:** **Unconditionally Stable** ($\forall \sigma$). Time step size is governed by accuracy rather than stability limits.

### Boundary Conditions

* **Inflow ($j=0$):** Dirichlet ($U_0^{n+1} = 0$)
* **Outflow ($j=N_x-1$):** Implemented directly into row $N_x-2$ of matrix $A$ as $U_{N_x-1}^{n+1} = U_{N_x-2}^{n+1}$, yielding $(1 + \frac{\sigma}{2}) U_{N_x-2}^{n+1} - \frac{\sigma}{2} U_{N_x-3}^{n+1} = U_{N_x-2}^n$.

---

## Summary Comparison Table

| Method | Order (Time, Space) | Type | Stability Limit | Dissipation | Dispersion |
| --- | --- | --- | --- | --- | --- |
| **FTCS** | $(\mathcal{O}(\Delta t), \mathcal{O}(\Delta x^2))$ | Explicit | Unstable ($\forall \sigma$) | Negative (Blows up) | Severe |
| **Upwind** | $(\mathcal{O}(\Delta t), \mathcal{O}(\Delta x))$ | Explicit | $\sigma \le 1$ | High (Damping) | Moderate |
| **Lax-Friedrichs** | $(\mathcal{O}(\Delta t), \mathcal{O}(\Delta x^2/\Delta t))$ | Explicit | $\sigma \le 1$ | Very High | High |
| **Lax-Wendroff** | $(\mathcal{O}(\Delta t^2), \mathcal{O}(\Delta x^2))$ | Explicit | $\sigma \le 1$ | Low | High (Trailing oscillations) |
| **Leapfrog** | $(\mathcal{O}(\Delta t^2), \mathcal{O}(\Delta x^2))$ | Explicit | $\sigma \le 1$ | None ($\vert{}G\vert{}=1$) | High (Trailing oscillations + Decoupling) |
| **BTCS** | $(\mathcal{O}(\Delta t), \mathcal{O}(\Delta x^2))$ | Implicit | Unconditionally Stable | Moderate to High | High at large $\sigma$ |
<img width="1000" height="600" alt="Figure_1 2" src="https://github.com/user-attachments/assets/7ec35293-4e0e-4aea-bc67-09b90068f3b1" />
