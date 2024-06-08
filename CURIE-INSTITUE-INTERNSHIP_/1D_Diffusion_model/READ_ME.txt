In this folder you will find the different versions of my codes where the different diffusion models of hormones are.

In Diffusion_model_1D_v1 we plot 
- the source term shape in space
- the initial conditions at boundaries for the hormone concentration
- the concentration of hormones in the middle of the worm in time
- the concentration of hormones at stationnary state computed with the given boundary conditions in space
- the concentration of hormones at stationnary state of a gaussian analytical model where the worm would be an infinite worm with the given boundary conditions in space
- the flux in space at stationnary state
- the flux in time in the middle of the worm
- the flux in time at the edges

In Diffusion_model_1D_v2 we plot:
- the source term shape in space
- the initial conditions at boundaries for the hormone concentration
- the concentration of hormones in the middle of the worm in time
- the concentration of hormones at stationnary state computed with the given boundary conditions in space
- the concentration of hormones at stationnary state of a gaussian analytical model where the worm would be an infinite worm with the given boundary conditions in space
- the standard deviation of the concentration of hormones curve for every timepoint => sigma(t) = f([Hormones](x;t))
- the flux in space in all the worm at the end of the simulation
- the flux in time in the middle of the worm
- the flux in time at the edges

In Diffusion_model_1D_v2_line_cells we plot almost the same thing as the previous one but it was an intermediate version where errors where about to be solved but wern't at the end.
This verrsion isn't working but the core of the code is the same.

In Diffusion_model_1D_v3 we plot:
- the source term shape in space
- the initial concentration of hormones in space
- the concentration in the center in time
- the evolution of the concentration at every timepoint - Numeric
- the evolution of the concentration at every timepoint - Analytic - Gaussian shape
- the standard deviation in time
- the flux in space at stationary state

In Diffusion_model_v3_without_last_modifs we plot:
- the source term shape in space
- the concentration at stationary state - Numeric

In Diffusion_model_1D_v4_Ncells_Shape_recognition we plot:
- the source term shape in space
- the initial concentration of hormones in space
- the concentration of hormones in the middle of the worm in time
- the evolution of the concentration at every timepoint - Numeric
- the evolution of the concentration at every timepoint - Analytic - Gaussian shape
- the standard deviation in time
- the flux in space at stationary state
- the error of the concentration distribution at every timepoinbt compared to an x² theoretical concentration distribution in space (this is computed according to the least square method (LSM))

In Diffusion_model_1D_v5_Ncells_nucleus_taken_into_account we plot:
- the source term shape in space
- the initial concentration of hormones in space
- the concentration of hormones in the middle of the worm in time
- the evolution of the concentration at every timepoint - Numeric
- the standard deviation in time
- the flux in space at stationary state
- the error of the concentration distribution at every timepoinbt compared to an x² theoretical concentration distribution in space (this is computed according to the least square method (LSM))


For the other three codes that are there, it is just an extract of all the previous one and they are more testing codes but still can be used so I let them here.