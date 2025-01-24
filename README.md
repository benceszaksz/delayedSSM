# delayedSSM: MATLAB code for the computation of invariant manifolds in time delay systems

This is a MATLAB package for the calculation of the spectral submanifold (SSM) and the corresponding reduced dynamics in time delay systems. 

The corresponding theoretical background is described in the following articles:

**[1] Szaksz, B., Orosz, G., Stepan, G. Spectral submanifolds in time delay Systems. Nonlinear Dyn. (2025) https://doi.org/10.1007/s11071-025-10902-0**

**[2] Szaksz, B., Orosz, G., Stepan, G. Reduction to spectral submanifolds in guided car-following with time delay. IFAC-PapersOnLine 58(27): 67-72 (2024) https://doi.org/10.1016/j.ifacol.2024.10.301**

The package can be installed by running the install.m file, which adds the "function" folder to the MATLAB path. Similarly, the uninstall.m file removes the function folder from the MATLAB path.

Currently, there are 3 demos within the demos folder. Each case study is in a separate folder, which contains the corresponding equation of motion in the file "EoM_..." and the main file(s) called "SSM_main_...".

For the sake of simplicity, separate files are created for the case of a real dominant eigenvalue and for the case of a pair a complex conjugate dominant eigenvalues. The toolbox can handle these cases in one file, only the visualization should be separated accordingly.

Main features:
- red_dyn_style
One can specify the style of the reduced dynamics. It can be either 'graph', 'normal_form' (only in case of complex conjugate roots), or 'manual'

- nlin_type
One can specify the type of nonlinearity to make the calculations faster. It can be either 'non-delayed', 'delayed' or 'combined.

- SSMorder
Order of the SSM and the corresponding reduced dynamics. The toolbox can handle high order cases, however, because of the increasing computational time and decreasing numerical accuracy, it is not suggested to go above O(7).

The current demos are:
- Scalar: 
Nonlinear scalar delay differential equation with symmetric nonlinearity. A limit cycle exists for appropriate parameters, which is obtained from the reduced dynamics if one runs the SSM_main_scalar_complex.m file.

- Scalar_bistab:
Nonlinear scalar delay differential equation with nonsymmetric nonlinearity. The bifurcation diagram has a fold point which can be predicted. It is worth running the code with SSMorder=7. 

- Guided_control:
Nonlinear guidance of human driver via an automated vehicle. It yields a 3-dimensional delay differential equation with both second and third-order nonlinearities.
