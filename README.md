# LSC MAIN PROJECT: RESISTIVE MHD SOLVER
Built a 2D resistive MHD solver in C++ and investigated a range of standard tests that demonstrate magnetic reconnection

## Code Details
- Finite Volume Scheme using a HLLD (5 wave approximate Riemann) conservative flux update scheme
- Uses an additional generalised lagrange multipler (GLM) wave to enforce divergence clearing (divB=0)
- Performs explicit RK2 updates for resistive source terms with optional subcycling


## Orszag-Tang Vortex Test
- Small but finite resistivity (η) -> magnetic reconnection at domain centre
- Visible as formation and ejection of plasmoids
- At larger resistivities increased diffusivity suppresses this

<table>
  <tr>
    <td>
      <img src="ReadMeFigs/eta_10minus4.png" width="300"/><br/>
      <p align="center"><em>Figure 1: η=10^(-4) no plasmoid formation visible </em></p>
    </td>
    <td>
      <img src="ReadMeFigs/eta_10minus5.png" alt="Orszang-Tang Animation" width="300"/><br/>
      <p align="center"><em>Figure 2: η=10^(-4) magnetic reconnection at centre visible</em></p>
    </td>
  </tr>
</table>




## GEM Test


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


## Double Current Sheet


<td>
<img src="FIGURES/solution.png" alt="Solov'Ev Plot" width="300"/>
<p align="left"><em>Figure 3: Sparse Linear solver solution to Solov'Ev equation</em></p>
</td> 





  

