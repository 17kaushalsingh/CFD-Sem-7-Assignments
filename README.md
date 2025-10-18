# 2D Streamfunction-Vorticity CFD Solver
## Technical Documentation

---

## 1. Problem Statement

### 1.1 Physical Problem
This solver simulates incompressible, viscous fluid flow in a 2D rectangular channel domain. The physical configuration consists of:

- **Domain dimensions**: Length L = 10 m, Height H = 1 m
- **Flow characteristics**: Inlet velocity U₀ = 1 m/s, Reynolds number Re = 1000
- **Flow type**: Laminar channel flow with parabolic inlet profile

The problem represents a classic computational fluid dynamics test case used to validate numerical methods for incompressible Navier-Stokes equations.

### 1.2 Governing Equations
The incompressible Navier-Stokes equations in 2D are reformulated using the streamfunction-vorticity (ψ-ω) approach:

**Vorticity Transport Equation:**
```
∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = (1/Re)(∂²ω/∂x² + ∂²ω/∂y²)
```

**Streamfunction-Vorticity Relation:**
```
∇²ψ = ω
```
or equivalently:
```
∂²ψ/∂x² + ∂²ψ/∂y² = ω
```

**Velocity-Streamfunction Relations:**
```
u = ∂ψ/∂y
v = -∂ψ/∂x
```

where:
- ω = vorticity (∂v/∂x - ∂u/∂y)
- ψ = streamfunction
- u, v = velocity components in x and y directions
- Re = Reynolds number (ρUL/μ)

### 1.3 Advantages of Streamfunction-Vorticity Formulation
This formulation offers several benefits:
- **Automatic continuity satisfaction**: The continuity equation (∇·u = 0) is automatically satisfied by the velocity-streamfunction relationship
- **No pressure term**: Eliminates the pressure field from the momentum equations, reducing computational complexity
- **Natural boundary conditions**: Wall conditions are naturally expressed in terms of streamfunction values

---

## 2. Numerical Method

### 2.1 Discretization Approach
The solver employs a **finite difference method** with:
- **Spatial discretization**: Second-order central differences
- **Temporal discretization**: Explicit forward Euler for vorticity transport
- **Grid**: Uniform rectangular mesh with nx × ny points

### 2.2 Grid Generation
The computational domain is discretized uniformly:
- x-direction: nx = 201 points, spacing dx = L/(nx-1)
- y-direction: ny = 51 points, spacing dy = H/(ny-1)
```
x ∈ [0, L], with dx = 0.05 m
y ∈ [0, H], with dy = 0.02 m
```

### 2.3 Solution Algorithm
The solver uses a **segregated approach** with operator splitting:

**Time-stepping loop:**
1. **Solve Poisson equation for ψ**: Given ω, solve ∇²ψ = ω iteratively
2. **Compute velocities**: Calculate u and v from ψ using finite differences
3. **Advance vorticity**: Update ω using the vorticity transport equation
4. **Apply boundary conditions**: Update boundary values for both ψ and ω
5. **Check convergence**: Monitor residuals and iterate

---

## 3. Solution Methods

### 3.1 Streamfunction Solver (Poisson Equation)
The Poisson equation ∇²ψ = ω is solved using **Successive Over-Relaxation (SOR)** with point-wise Gauss-Seidel iteration.

**Discretized equation:**
```
ψᵢⱼ = [dx²(ψᵢ,ⱼ₊₁ + ψᵢ,ⱼ₋₁) + dy²(ψᵢ₊₁,ⱼ + ψᵢ₋₁,ⱼ) + dx²dy²ωᵢⱼ] / [2(dx² + dy²)]
```

**Algorithm:**
- Iterative method with relaxation parameter ω_relax = 1.0 (Gauss-Seidel)
- Convergence tolerance: 10⁻⁸
- Maximum iterations: 3000
- Updates are performed in-place for faster convergence

**Why this method?**
- Simple to implement and parallelize
- Robust convergence for elliptic problems
- Well-suited for rectangular domains
- Low memory requirements

### 3.2 Vorticity Advection-Diffusion Solver
The vorticity transport equation is advanced using **explicit Euler time integration** with central differences for spatial derivatives.

**Discretization:**

**Advection term:**
```
ADV = u·∂ω/∂x + v·∂ω/∂y
    ≈ uᵢⱼ(ωᵢ₊₁,ⱼ - ωᵢ₋₁,ⱼ)/(2dx) + vᵢⱼ(ωᵢ,ⱼ₊₁ - ωᵢ,ⱼ₋₁)/(2dy)
```

**Diffusion term:**
```
DIFF = (1/Re)∇²ω
     ≈ (1/Re)[(ωᵢ₊₁,ⱼ - 2ωᵢⱼ + ωᵢ₋₁,ⱼ)/dx² + (ωᵢ,ⱼ₊₁ - 2ωᵢⱼ + ωᵢ,ⱼ₋₁)/dy²]
```

**Time advancement:**
```
ωⁿ⁺¹ᵢⱼ = ωⁿᵢⱼ + Δt(-ADV + DIFF)
```

**Why explicit method?**
- Straightforward implementation
- No need for matrix solvers
- Suitable for moderate CFL numbers
- Time step automatically adjusted for stability

---

## 4. Stability Analysis

### 4.1 CFL Condition (Advection)
For the advection term, the **Courant-Friedrichs-Lewy (CFL)** condition must be satisfied:
```
CFL = |u|Δt/Δx + |v|Δt/Δy ≤ C_max
```

where C_max ≈ 1 for stability.

**Implemented constraint:**
```
Δt_conv = 0.5 × min(Δx, Δy) / max(u₀, 10⁻⁸)
```

This ensures CFL ≈ 0.5, providing a safety margin.

### 4.2 Diffusion Stability
For the diffusion term with explicit Euler, the **von Neumann stability** condition requires:
```
D_max = (Δt/Re) × (1/Δx² + 1/Δy²) ≤ 1/2
```

**Implemented constraint:**
```
Δt_diff = 0.25 × min(Δx², Δy²) × Re
```

This ensures the diffusion number D ≤ 0.25, providing stability with margin.

### 4.3 Final Time Step Selection
The code automatically selects:
```
Δt = min(Δt_conv, Δt_diff)
```

**For default parameters:**
- Δt_conv = 0.5 × 0.02 / 1.0 = 0.01 s
- Δt_diff = 0.25 × 0.0004 × 1000 = 0.1 s
- **Selected Δt = 0.01 s** (limited by convection)

---

## 5. Initial Conditions

### 5.1 Streamfunction Initialization
The streamfunction is initialized with a **fully developed parabolic profile** to accelerate convergence:

**Inlet velocity profile:**
```
u_inlet(y) = 6U₀(y/H)(1 - y/H)
```

This parabolic profile satisfies:
- No-slip at walls: u(0) = u(H) = 0
- Mean velocity: ∫u dy / H = U₀

**Streamfunction from velocity:**
```
ψ(y) = ∫₀ʸ u_inlet(η) dη
```

Computed using trapezoidal integration:
```
ψⱼ = ψⱼ₋₁ + 0.5(uⱼ + uⱼ₋₁)(yⱼ - yⱼ₋₁)
```

**Applied everywhere:**
The inlet profile is extended throughout the domain as the initial guess, which represents a reasonable starting point.

### 5.2 Vorticity Initialization
Initial vorticity is set to **zero in the interior** and then boundary conditions are applied immediately before the first iteration.

---

## 6. Boundary Conditions

### 6.1 Streamfunction Boundary Conditions

**Bottom and Top Walls (y = 0, H):**
```
ψ(x, 0) = 0
ψ(x, H) = ψ_inlet(H) = constant
```

These ensure no flow through walls (u·n = 0).

**Inlet (x = 0):**
```
ψ(0, y) = ψ_inlet(y)
```

Prescribed parabolic profile.

**Outlet (x = L):**
Streamfunction naturally advects out; typically:
```
∂ψ/∂x = 0  (Neumann condition)
```

### 6.2 Vorticity Boundary Conditions

**Bottom and Top Walls (y = 0, H):**
Wall vorticity enforces no-slip condition using **Thom's formula**:
```
ω(x, 0) = -2[ψ(x, 1) - ψ(x, 0)] / dy²
ω(x, H) = -2[ψ(x, ny-2) - ψ(x, H)] / dy²
```

This ensures the velocity gradient at the wall matches the no-slip requirement.

**Inlet (x = 0):**
Vorticity consistent with parabolic inlet profile:
```
ω_inlet(y) = -du_inlet/dy
```

For the parabolic profile:
```
ω_inlet(y) = -6U₀/H[1 - 2y/H]
```

Computed using `np.gradient()` for accuracy.

**Outlet (x = L):**
**Convective boundary condition** to allow vortices to exit smoothly:
```
∂ω/∂t + U_local·∂ω/∂x = 0
```

Discretized using upwind differences:
```
ω(L, y)ⁿ⁺¹ = ω(L-dx, y)ⁿ - (U_local·Δt/Δx)[ω(L-dx, y) - ω(L-2dx, y)]
```

where U_local = u(L-dx, y) is the local velocity.

**Corner Points:**
Averaged from adjacent boundaries:
```
ω(0, 0) = 0.5[ω(1, 0) + ω(0, 1)]
```
and similarly for other corners.

---

## 7. Convergence Criteria

### 7.1 Inner Iteration (Poisson Solver)
The streamfunction solver iterates until:
```
max|ψⁿ⁺¹ - ψⁿ| < tol_ψ = 10⁻⁸
```

### 7.2 Outer Iteration (Time Stepping)
The simulation continues until vorticity field converges:
```
||ωⁿ⁺¹ - ωⁿ||₂ / (||ωⁿ||₂ + ε) < tol_ω = 10⁻⁶
```

where ε = 10⁻¹⁴ prevents division by zero.

This indicates the solution has reached steady-state.

---

## 8. Output and Post-Processing

### 8.1 Solution Files
Output is written in **Tecplot ASCII format** at specified intervals:
- Files: `solution_001.dat`, `solution_002.dat`, ...
- Contains: x, y, ψ, ω, u, v at all grid points
- Frequency: Every 10 iterations (default)

### 8.2 Residual Logging
Convergence history is logged to `omegarsd.dat`:
```
Iter  Time[s]  OmegaResidual
   1   0.0100  1.234567e-02
   2   0.0200  9.876543e-03
  ...
```

### 8.3 Final Output
Upon convergence or reaching maximum iterations:
- `solution_final.dat`: Final converged solution
- Console output: Convergence message and timing statistics

---

## 9. Computational Performance

### 9.1 Complexity
- **Per iteration**: O(nx × ny × n_inner)
- **Total**: O(n_outer × nx × ny × n_inner)

For default parameters:
- Grid: 201 × 51 = 10,251 points
- Typical iterations: ~100-500 outer, ~100-500 inner per outer
- Runtime: Minutes on modern CPU

### 9.2 Optimization Opportunities
Current implementation priorities:
1. **Code clarity** over raw performance
2. **Stability** over speed
3. **Extensibility** for future enhancements

Potential improvements:
- Vectorization using NumPy operations
- Multigrid methods for Poisson solver
- Implicit time stepping for larger time steps
- Parallel processing (OpenMP/MPI)

---

## 10. Validation and Testing

### 10.1 Expected Results
For the channel flow configuration:
- Parabolic velocity profile should develop and maintain
- Vorticity concentrated near walls
- Steady-state reached after transient development

### 10.2 Key Dimensionless Parameters
- **Reynolds number (Re = 1000)**: Laminar flow regime
- **Aspect ratio (L/H = 10)**: Long channel allows development
- **Grid resolution**: Adequate for Re = 1000

### 10.3 Physical Consistency Checks
- Mass conservation: ∫u dy = constant
- Maximum velocity: u_max ≈ 1.5 U₀
- Wall vorticity: Proportional to wall shear stress

---

## 11. Usage Example
```bash
# Run with default parameters
python solver.py

# Custom grid and Reynolds number
python solver.py --nx 301 --ny 101 --Re 500

# Finer time step
python solver.py --dt 0.005 --itmax 5000

# More frequent output
python solver.py --outfreq 5
```

---

## 12. References and Further Reading

1. **Streamfunction-Vorticity Methods**: Classical approach for 2D incompressible flow
2. **Finite Difference Methods for PDEs**: Theoretical foundation
3. **Computational Fluid Dynamics**: Numerical techniques and stability analysis
4. **Thom's Formula**: Boundary vorticity calculation for no-slip walls

---

## Appendix: Key Parameters Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| nx | 201 | Grid points in x |
| ny | 51 | Grid points in y |
| L | 10.0 m | Domain length |
| H | 1.0 m | Domain height |
| Re | 1000 | Reynolds number |
| U₀ | 1.0 m/s | Inlet velocity scale |
| Δt | auto | Time step (stability-limited) |
| tol_ψ | 10⁻⁸ | Poisson solver tolerance |
| tol_ω | 10⁻⁶ | Convergence tolerance |
| itmax | 2000 | Maximum iterations |

---

**Document Version**: 1.0  
**Date**: October 2025  
**Software**: Python 3.x with NumPy