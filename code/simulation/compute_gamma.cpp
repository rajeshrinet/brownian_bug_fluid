//23/01/2020 CP: This code uses different Utot, without diffusion and birth/death, and outputs the position of pairs of particles. This should allow us to compute the gamma parameter to solve eq. (2)

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "basic_particle.h"
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

gsl_rng *rgslbis2 = gsl_rng_alloc(gsl_rng_mt19937);

//Define constant
extern const double pi=3.14159265;
extern const double Delta=0.0; //diffusion
extern const double Lmax=1.0; //size of the grid
extern const double k=2*pi/Lmax; 
extern const int size_pop=800; //Number of pairs of particles: the total size of the population is 2*size_pop 
extern const int tmax=100; //length of the simulation. The system stabilizes very quickly, especially for high values of Utot, we do not need more than 100 points
extern const double proba_death=0.5; //Death and birth probability
extern const double proba_repro=0.5;

using namespace std;

//Basic functions to compute the distance between two particles
double distance(basic_particle p1, basic_particle p2)
{
	double x1,x2,y1,y2,dist;
	x1=p1.get_x();
	x2=p2.get_x();
	y1=p1.get_y();
	y2=p2.get_y();
	dist=pow(x1-x2,2.0)+pow(y1-y2,2.0);
	return dist;
}

int main()
{
	int i,j,t,tmp_t;
	double a_x,a_y,phi,theta,a_n,xi,dxi,pow_min,pow_max,dpow,pow_i,pcf;
	std::vector<basic_particle> Part_table_parent,Part_table_children; //Here, there are two populations: each parent particle has a counterpart/children, from which it has been separated by 10^(-7) at the beginning of the simulation. This allows us to identify easily the pairs between which we want to compute distance
	std::vector<double> Utot_list{ 0.0, 0.1, 0.5,2.5 };
	std::ofstream f0,f1;

	//Open the file in which we will have the x, y, parent of each particle
	f0.open("Spatial_distribution_particle_compute_gamma.txt");
	f0.precision(8);
	for (double Utot : Utot_list) 
	{
	//Initialize
	for(i=0; i < size_pop; i++)
	{
		a_x=gsl_rng_uniform(rgslbis2);
		a_y=gsl_rng_uniform(rgslbis2);
		Part_table_parent.push_back(basic_particle(a_x,a_y,a_y,i));
                f0<< Utot <<";";
                f0<< 0 <<";";
                f0<< Part_table_parent[i].get_x()  << ';';
                f0<< Part_table_parent[i].get_y()  << ';';
                f0<< Part_table_parent[i].get_yfirst()  << ";";
                f0<< "P"  << ";";
                f0<< Part_table_parent[i].get_firstparent()  << std::endl;
		Part_table_children.push_back(basic_particle(a_x,a_y+pow(10,-7),a_y+pow(10,-7),i)); //The particles should be separated by 10^-7: I only put this distance in the y axis.
                f0<< Utot <<";";
                f0<< 0 <<";";
                f0<< Part_table_children[i].get_x()  << ';'; 
                f0<< Part_table_children[i].get_y()  << ';';
                f0<< Part_table_children[i].get_yfirst()  << ";";
                f0<< "C"  << ";";
                f0<< Part_table_children[i].get_firstparent()  << std::endl; 
	}

	//Run the simulation
	for(t=0;t<=tmax;t++)
	{
	printf("TIME=%d, SIZE=%u\n",t,Part_table_parent.size());
	//Compute the phase in x and y for the turbulent flow from Pierrehumbert. These phases are common to each particle as they correspond to a unique flow
	phi=gsl_rng_uniform(rgslbis2)*2*pi;
	theta=gsl_rng_uniform(rgslbis2)*2*pi;

	//For each particle, diffusion, then advection
	for(j=0; j<Part_table_parent.size(); j++)
	{
		Part_table_parent[j].pierrehumbert_flow(Utot, k, phi,theta, Lmax);
		Part_table_children[j].pierrehumbert_flow(Utot, k, phi,theta, Lmax);
                
		f0<< Utot <<";";
                f0<< t <<";";
                f0<< Part_table_parent[j].get_x()  << ';';
                f0<< Part_table_parent[j].get_y()  << ';';
                f0<< Part_table_parent[j].get_yfirst()  << ";";
		f0<< "P"  << ";";
                f0<< Part_table_parent[j].get_firstparent()  << std::endl;
                f0<< Utot <<";";
                f0<< t <<";";
                f0<< Part_table_children[j].get_x()  << ';';
                f0<< Part_table_children[j].get_y()  << ';';
                f0<< Part_table_children[j].get_yfirst()  << ";";
                f0<< "C"  << ";";
                f0<< Part_table_children[j].get_firstparent()  << std::endl;


	} //end j, i.e. the particle mvt
	} //end t, i.e the whole simulation
	 Part_table_children= std::vector<basic_particle>(); //Reinitialize the particles
	 Part_table_parent= std::vector<basic_particle>();
	}
	//End of the simulation on U
	f0.close();
	return 0;
}

int main();
