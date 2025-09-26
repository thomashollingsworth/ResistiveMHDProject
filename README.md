# LSC MAIN PROJECT: RESISTIVE MHD SOLVER
Built a 2D resistive MHD solver in C++ and investigated a range of standard tests that demonstrate magnetic reconnection

## Code Details
- Finite Volume Scheme using a HLLD (5 wave approximate Riemann) conservative flux update scheme
- Uses a generalised lagrange multipler (GLM) wave to enforce divergence clearing (divB=0)
- Performs RK2 updates for resistive source terms with optional subcycling


## Solving Euler equations in 1D & 2D
  - Using Lax Friedrichs + Richtymer fluxes and First ORder CEntered (FORCE) scheme
  - Incorporated SLIC (slope limiting) with MinBee limiter

## MHD Solvers in 1D & 2D
  - Created Godunov solvers using HLL and HLLC (2 wave and 3 wave approximations of exact Riemann problem)
  - Incorporated Van-Leer slope limiting -> MUSCL-Hancock scheme
  - Implemented Divergence Cleaning
  - Validation of sovlers using Brio-Wu (shock tube) tests in 1 and 2D
  - Further validation using Orszang Tang and Kelvin-Helmholtz tests

<table>
  <tr>
    <td>
      <img src="FIGURES/KHAni.gif" width="300"/><br/>
      <p align="center"><em>Figure 1: Animation of sqrt(Bx^2 + By^2) / Bz (poloidal to toroidal field strength ratio) for Kelvin-Helmholtz Instability test using MUSCL-Hancock Scheme</em></p>
    </td>
    <td>
      <img src="FIGURES/HighResDensityAni.gif" alt="Orszang-Tang Animation" width="300"/><br/>
      <p align="center"><em>Figure 2: Animation of Density for Orszang-Tang vortex using MUSCL-Hancock Scheme</em></p>
    </td>
  </tr>
</table>


## Linear Solver for Solov'Ev equation
  - Using C++ eigen library to solve sparse linear Solov'ev equation in 2D (axisymmetric case)

<td>
<img src="FIGURES/solution.png" alt="Solov'Ev Plot" width="300"/>
<p align="left"><em>Figure 3: Sparse Linear solver solution to Solov'Ev equation</em></p>
</td> 





  

