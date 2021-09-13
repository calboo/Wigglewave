# Wigglewave

WiggleWave is a FORTRAN code that uses an RK4 finite difference method to solve the linearised governing equations for a torsion Alfv&egrave;n wave propagating in a plasma with negligible plasma beta and in a force-free axisymmetric magnetic field with no azimuthal component embedded in a high density divergent tube structure. The solutions calculated are the perturbations to the velocity, v and to the magnetic field, b. All variables are calculated over a uniform grid in radius r and height z.

## Usage

TAWAS can be run simply using IDL. 

The code requires no inputs but the user is able to change the problem parameters at the top of the code, the parameters are listed in the first table below. Under the parameters the user can set the variable 'plotting' which decides which plots are shown by the script. The options for plotting are given in the code. The user can chose the location to save the outputs by changing the variables save_dir and run_name listed just below the parameters.

The code outputs include solutions for the velocity perturbation, the magnetic field perturbation and the wave energy flux . The outputs are saved as .sav files and are listed in the second table below.

## Input Parameters

The parameters that can be changed are at the begining of the code. These parameters are:

| Parameter | Description |
| --- | --- |
| b0     | background magnetic field strength                               |
| rho0   | characteristic density                                           |
| rsize  | number of cells in redial direction                              |
| zsize  | number of cells in vertical direction                            |
| rdim   | dimension of domain in the radial direction                      |
| zdim   | dimension of domain in the vertical direction                    |
| H      | magnetic scale height                                            |
| alpha  | &alpha; parameter, defines density scale height through 	&alpha; = H/H<sub>&rho;</sub>    |
| r0     | radius of central higher density tube                            |
| zeta   | density contrast between tube centre and background density      |
| u0     | Alfv&egrave;n wave velocity amplitude in ms<sup>&-1;</sup>       |
| period | wave period in seconds                                           |
| t0     | wave phase of solution as time in seconds (usually left as zero) |

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
