# Wigglewave

WiggleWave is a FORTRAN code that uses an RK4 finite difference method to solve the linearised governing equations for a torsion Alfv&egrave;n wave propagating in a plasma with negligible plasma beta and in a force-free axisymmetric magnetic field with no azimuthal component embedded in a high density divergent tube structure. The solutions calculated are the perturbations to the velocity, v and to the magnetic field, b. All variables are calculated over a uniform grid in radius r and height z. An in-depth explanation of the code can be found in the accompanying PDF document Wigglewave equations.

## Usage

WiggleWave can be run simply compiled and run using gfortran. 

The code requires no input files but the user is able to change the problem parameters in the Constants module, the parameters are listed in the first table below.

The code outputs include solutions for the velocity perturbation, the magnetic field perturbation and wave envelopes for these perturbations. The outputs are saved as .dat files and are listed in the second table below.

To calculate the total wave energy flux from the wave envelopes the user must first convert the output files from .dat to .sav files using the IDL script makesavs.pro and then use the IDL script wave_enrgy.pro on the .sav files generated.

## Input Parameters

The parameters that can be changed are at the begining of the code. These parameters are:

| Parameter | Description |
| --- | --- |
| B0        | background magnetic field strength                               |
| rho0      | characteristic density                                           |
| nr        | number of cells in redial direction                              |
| nz        | number of cells in vertical direction                            |
| rmin      | minimum radius of domain                                         |
| rmax      | maximum radius of domain                                         |
| zmin      | minimum height of domain                                         |
| zmax      |  maximum height of domain                                        |
| t_end     |  simulation run time                                             |
| t_interval|  time interval between outputs                                   |
| t0        |  rampup time                                                     |
| save_dir  |  directory to save outputs to                                    |
| H      | magnetic scale height                                               |
| visc   | kinematic viscosity                                                 |
| period | wave period in seconds                                              |
| alpha  | &alpha; parameter, defines density scale height through 	&alpha; = H/H<sub>&rho;</sub>    |
| zeta   | density contrast between tube centre and background density         |
| u0     | Alfv&egrave;n wave velocity amplitude in ms<sup>&-1;</sup>          |
| r0     | radius of central higher density tube                               |
| omega  | the Alfv&egrave;n  wave frequency, currently defined based on the period |
| topdamp  | logical for top boundary damping|
| outdamp  | logical for outer boundary damping |
| restart  | logical for whether to load from a restart file |
| v_in     | file names for the velocity restart file |
| b_in     | file names for the magnetic field restart file|
| restart_time | the simulation time for the restart files |
| last_output  | the output index for the restart files |

## Outputs

| Output | Format | Description |
| --- | ----- | --- |
| br           | 2D array              | background magnetic field in radial direction        |
| bz           | 2D array              | background magnetic field in vertical direction      |
| rho          | 2D array              | density across domain                                |
| phi          | 2D array              | curvilinear coorinate along field lines &phi;        |
| psi          | 2D array              | curvilinear coorinate along field lines &psi;        |
| hparam       | 2D array              | phase parameter h , see below for more details       |
| W            | 2D complex array      | wave parameter W , see below for more details        |
| v_env        | 2D complex array      | envelope of the velocity perturbation                |
| b_env        | 2D complex array      | envelope of the magnetic field perturbation          |
| v            | 2D complex array      | azimuthal velocity perturbation                      |
| bwave        | 2D complex array      | azimuthal magnetic field perturbation                |
| sigma        | 2D array              | WKB parameter &sigma;, see below for more details    |
| zscale       | 1D array              | heights at which each magnetic surface intersects the z-axis, in metres       |
| hscale       | 1D array              | heights at which each magnetic surface intersects the z-axis, in units of H   |
| en_lvl       | 1D array              | wave energy flux across each magnetic surface in Watts                        |
| en_lvl_norm  | 1D array              | wave energy flux normalised by the wave energy flux at the lowest surface     |

## makesavs.pro

## wave_energy.pro
