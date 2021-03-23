## Description

* `find_gamma.r` uses results from the script `compute_gamma.cpp` to compute the stretching parameter ![\gamma](https://latex.codecogs.com/svg.latex?\gamma), the timescale for particle separation distance, under different advection values.  
* `visualisation_Fig*.r` use the results from `Spatial_distribution_particle*.txt` to produce Fig. 1 and 2 of the original paper, i.e. the spatial repartition of the bugs with and without advection and/or birth/death. They output the figures `spatial_distribution_Fig*.png`
*  `pcf_Fig3.r` uses the `pcf_particle*.txt` files to compare the simulated partial correlation function values to the analytical ones. Analytical formula are justified in the folder `article`  
