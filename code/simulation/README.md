## Description

* `basic_particle` files define the class describing an individual organism. 
* `main_Fig*.cpp` files reproduce the simulations described in the original paper (numbers of the figure remain the same in the original paper, the replication, and the name of the files). They differ from one another mainly by the number of particles, value of diffusion, or phenomena taken into account. 
* `compute_gamma.cpp` is used to approximate the value of $\gamma$ by computing the separation between particles as a function of time (details are given in the Methods of the original paper). 

Executable files can be produced with the command `make` once `makefile` has been adapted.

Extraction of relevant values (stretching parameter or pair densitiy) and figures showing the results of these simulations can be found in the folder `../figure`.
