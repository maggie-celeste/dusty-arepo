# Arepo public version
====================

AREPO is a massively parallel code for gravitational n-body 
systems and hydrodynamics, both on Newtonian as well as 
cosmological background. It is a flexible code that can be 
applied to a variety of different types of simulations, offering 
a number of sophisticated simulation algorithms. An description 
of the numerical algorithms employed by the code is given in the 
original code papers (Springel 2010, MNRAS, 401, 791; 
Pakmor et al. 2011, MNRAS, 418, 1392; Pakmor and Springel 2013, 
MNRAS, 432, 176; Pakmor et al. 2016, MNRAS,455,1134) and the 
release paper of this version (Weinberger et al. 2019). 

A user guide can be found under `/documentation`, which also 
includes a 'getting started' section, which is recommended for 
new users. An html version of the user guide can be created using
sphinx (https://www.sphinx-doc.org) by typing

    cd ./documentation/
    make html
    
and displayed by opening `./documentation/build/html/index.html`.

A full version of the user guide is also available on the Arepo 
homepage.

## dusty-arepo

This branch of AREPO has been modified to include dust, modified Lombardi cooling, and a torque-free-sink particle. These changes make AREPO well-suited for protoplanetary disc simulations, but may also be useful in other cases. The additional compile time options needed to include dust, modified Lombardi cooling, and/or the torque-free-sink are outlined below.

## DUST

This code includes an algorithm for dust as a cold, pressureless fluid (Celeste, M et al, submitted). The initial conditions file must include dust masses ("DustMasses") and velocities ("DustVelocities") for each particle in the particle 0 (gas) dataset -- see the create.py file in examples/dustywave for an example of how to set this up.

The following compile-time options may be added to the Config.sh file:

> #DUST_INCLUDE

Basic flag for dust. Without this option, dust will not be included in the simulation.

> #DUST_STOKES = ST

An option for drag that sets the Stokes number, $St= \Omega \cdot t_s$, to be constant everywhere, equal to ST. Here $\Omega$ is the local angular frequency and $t_s$ is the single particle stopping time.

> #DUST_K = X

An option for drag that sets the drag coefficient, $K$, equal to X. $K$ is defined as:

$$K = \frac{1}{t_s (\rho_g + \alpha \rho_d)}$$

where $\rho_g$ is the gas density, $\rho_d$ is the dust density, and $\alpha$ is a flag for feedback of the dust on the gas.

> #DUST_SIZE = S

An option for drag in the Epstein regime. The dust grain size, $s$ is set equal to S. In this regime, the stopping time is then given by:

$$t_s = \frac{\rho_{grain} s_{grain}}{\rho_g c_g}$$

where $\rho_{grain}$ is the density of a dust grain and $c_g$ is the thermal speed of the gas.

> #DUST_RHO_GRAIN = Y

Option to set the dust grain density for the Epstein drag regime (doesn't do anything otherwise). Defaults to $2 gcm^{-3}$ if not set, or Y $g cm^{-3}$ if set.

> #DUST_FB

A flag to turn on or off feedback of dust onto gas. Not set => no feedback.

## COOLING

There are two cooling regimes available in this branch (along with the primordial cooling option of the master branch). Both options are available for a disc orbiting a single star, which must be included in the initial conditions. (See the general AREPO guide for how to initialise star particles).

# Modified Lombardi Cooling

Mdified Lombardi cooling is described in detail Young, A.K., Celeste, M., et al 2024, available here: https://academic.oup.com/mnras/article/531/1/1746/7671146?login=false 

To enable this version of cooling, use the following compile-time option:
> #MOD_LOMBARDI_COOLING

This is a flag to turn on modified Lombardi cooling. If this option is not present, the cooling will not be enabled.

> LombardiOpacityFile         ~/route/to/file/opacity_table.txt
> 
> LombardiPseudoOpacityFile   ~/route/to/file/pseudo_mean_opacity_table.txt

It is also necessary to input the location of the relevant opacity table files in the param.txt file as above, which is used to determine the opacity of a particle given a density and temperature. Two versions of these files have been included as .txt files in the src/cooling folder, with values generated using code available here: https://www.physics.unlv.edu/~zhzhu/Opacity.html -- one with ice included and one without. If you would like to generate your own opacity table, the file src/cooling/opacity_generate.cpp will be helpful.

> #OUTPUT_SCALEHEIGHT            

This is a flag that tells the code to output the scaleheight computed during the modified lombardi cooling routine. It is helpful if you want to calculate the cooling rate for particles in snapshots of a finished simulation.

# Beta Cooling

Beta cooling is a commonly used cooling regime in protoplanetary disc simulations. Proposed originally by Gammie (2001), the cooling timescale $t_c = \Beta / \Omega$. $\Beta$ may be set as a constant or allowed to vary with radius.

> #BETA_COOLING = B

Option to turn Beta-cooling on, where $\Beta$ = B.

> #BETA_POWER = P

Option to allow $\Beta$ to vary with radius from the central star as $\Beta = \Beta (R / R_0)^P$. P defaults to 0 if not set.

> #BETA_R0 = R

Option to allow $\Beta$ to vary with radius from the central star as $\Beta = \Beta (R / R_0)^P$. R defaults to 1 AU if not set. Make sure you set R with consistent units.

## Sink particle

To support simulations of protoplanetary discs, a torque-free-sink particle as described in Dempsey et al 2020 has been imlemented.

> #TORQUE_FREE_SINK

Turn on sink, centered on the first star particle. See the general AREPO guide for setting star particles.

The following options default to values that are appropriate for a protoplanetary disc of R_in 5AU, R_out 100 AU, with 2 million particles. Setting NPHI too low can result in the sink expanding and swallowing up the disc. Too high suffers a computational cost. The actual sink radius should be about half of the inner disc edge or less, to avoid swallowing the disc.

> #SINK_RADIUS

the radius inside of which sink particles are accreted onto the sink (not set=2AU)

> #INNER_RADIUS

the inner edge of the disc -- cells will gradually increase in size between inner_radius and sink_radius, but only be accreted inside sink_radius (see Dempsey et al 2021 for further info) (not set=5AU)

> #NPHI

How many cells across the circumference of the INNER EDGE? (not set=800)

> #ACCRETION_CONSTANT

the constant used to determine the sink rate, where rate(R=0) = const * omega (not set=10)

# Notes

When all of the above config options are disabled, this code is identical to the public AREPO master version (as of Sept 1st 2024).
The master repository of the AREPO code is available here:
https://gitlab.mpcdf.mpg.de/vrs/arepo/

Please report any bugs found, or get in touch if you have suggestions for improving the code or a feature you would like to see.