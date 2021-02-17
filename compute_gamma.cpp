//13/01/2021 CP: try and reproduce Young et al. (2001) code for the bug model
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
extern const int size_pop=800; //initial number of particles. Here, this the number of pairs. 
extern const int tmax=100; //length of the simulation. The system stabilizes very quickly, especially for high values of Utot, we do not need more than 100 points (even 50? Just checking longer for Utot=100)
extern const double proba_death=0.5; //Death and birth probability
extern const double proba_repro=0.5;

/*extern const unsigned seed=21; //random generator
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> distribution_unif_mv (0.0,Lmax);
std::uniform_real_distribution<double> distribution_unif_bio (0.0,1.0);
std::normal_distribution<double> distribution_normal (0.0,Delta);
std::uniform_real_distribution <double> dis_2pi(0.0,2*pi);
*/

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
	std::vector<basic_particle> Part_table_parent,Part_table_children;
	std::vector<double> Utot_list{ 0.0, 0.1, 0.5,2.5 };
	//std::vector<double> Utot_list{2.5 };
	std::ofstream f0,f1;

	dxi=0.001;
	pow_min=-1+log10(Delta);pow_max=5+log10(Delta); //These are the limits in Fig. 3 of Young et al. 2001
	dpow=0.25;

	//Open the file in which we will have the x, y, parent of each particle
	f0.open("Spatial_distribution_particle_compute_gamma.txt");
	f0.precision(8);
	for (double Utot : Utot_list) 
	{
	//Initialize
	for(i=0; i < size_pop; i++)
	{
		a_x=gsl_rng_uniform(rgslbis2);//distribution_unif_mv(generator);
		a_y=gsl_rng_uniform(rgslbis2);//distribution_unif_mv(generator);
		Part_table_parent.push_back(basic_particle(a_x,a_y,a_y,i));
                f0<< Utot <<";";
                f0<< 0 <<";";
                f0<< Part_table_parent[i].get_x()  << ';';
                f0<< Part_table_parent[i].get_y()  << ';';
                f0<< Part_table_parent[i].get_yfirst()  << ";";
                f0<< "P"  << ";";
                f0<< Part_table_parent[i].get_firstparent()  << std::endl;
		Part_table_children.push_back(basic_particle(a_x,a_y+pow(10,-7),a_y+pow(10,-7),i)); //Here, I am cheating. The particles should be separated by 10^-7: I do not distribute it between x and y.
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
	phi=gsl_rng_uniform(rgslbis2)*2*pi;//dis_2pi(generator);
	theta=gsl_rng_uniform(rgslbis2)*2*pi;//dis_2pi(generator);

	//For each particle, diffusion, then advection
	for(j=0; j<Part_table_parent.size(); j++)
	{
		Part_table_parent[j].pierrehumbert_flow(Utot, k, phi,theta, Lmax);
		Part_table_children[j].pierrehumbert_flow(Utot, k, phi,theta, Lmax);
//		if(((t>0)&(t%10==0))||(t==tmax))
//                {
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
//                }


	} //end j, i.e. the particle mvt
	} //end t, i.e the whole simulation
	 Part_table_children= std::vector<basic_particle>();
	 Part_table_parent= std::vector<basic_particle>();
	}
	//End of the simulation on U
	f0.close();
	return 0;
}

int main();
