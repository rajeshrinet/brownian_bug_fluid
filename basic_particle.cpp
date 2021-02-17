#include <stdio.h>
#include <stdlib.h>
#include "basic_particle.h"
#include <random>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

gsl_rng *rgslbis = gsl_rng_alloc(gsl_rng_mt19937);

//gsl_rng_set(rgslbis, 5);

using namespace std;

//Creator and destructor
basic_particle::basic_particle ()
{
  x=0;
  y=0;
  y_first=y;
}

basic_particle::~basic_particle()
{
}

basic_particle::basic_particle (double pos_x , double pos_y, double pos_yfirst, int first)
{
  x=pos_x;
  y=pos_y;
  y_first=pos_yfirst;
  first_parent=first;
}

//Use private variables
void basic_particle::set_x (double pos_x)
{
  x=pos_x;
}

void basic_particle::set_y (double pos_y)
{
  y=pos_y;
}

void basic_particle::set_yfirst (double pos_yfirst)
{
  y_first=pos_yfirst;
}

void basic_particle::set_firstparent (int first)
{
  first_parent=first;
}

double basic_particle::get_x()
{
  return x;
}

double basic_particle::get_y()
{
  return y;
}

double basic_particle::get_yfirst()
{
  return y_first;
}

int basic_particle::get_firstparent()
{
  return first_parent;
}

//Mvt of each particle
//Diffusion
void basic_particle::diffusion(double Delta, double Lmax)
{
//I do not know how to initialize the random function IN basic_particle.cpp and not only in main.cpp. Therefore, I need to change the seed, which is the current, varying time, with execution
	/*std::default_random_engine generator ((unsigned)time(NULL) );
	std::uniform_real_distribution<double> distribution_unif (0.0,Lmax);
	std::normal_distribution<double> distribution_normal (0.0,Delta);
*/


/*      
//This was the previous code, directly from Python. This does not seem right with the description in Young et al. Here, we have one random vector with modulus=1, then multiply its length by Delta. The total modulus is therefore the vaalue from the normal distribution
        double val_displacement,mv_x,mv_y,mov_x,mov_y;
	val_displacement=distribution_normal(generator);
        mv_x=distribution_unif(generator);
        mv_y=distribution_unif(generator);
        mov_x=mv_x/pow((pow(mv_x,2)+pow(mv_y,2)),0.5);
        mov_y=mv_y/pow((pow(mv_x,2)+pow(mv_y,2)),0.5);
        x=x+mov_x*val_displacement;
        y=y+mov_y*val_displacement;
*/

//Here, for each component of the vector, the diffusion adds a small increment with mean=0 and standard deviation=Delta. The total modulus of the increment is therefore dependent on two draws from the normal distribution, as opposed to what happened previously 
        double d_x,d_y;
//	d_x=distribution_normal(generator); 
	d_x=gsl_ran_gaussian(rgslbis,Delta);
//	d_y=distribution_normal(generator);
	d_y=gsl_ran_gaussian(rgslbis,Delta);
        x=x+d_x;
        y=y+d_y;
	check_boundaries(Lmax);
}

//Periodic conditions
void basic_particle::check_boundaries(double Lmax)
{
	x = fmod(x, Lmax);
	y = fmod(y, Lmax);
	if (x>Lmax)
	{
                x=x-Lmax;
	}
        else if (x<0)
	{
                x=x+Lmax;
	}

        if (y>Lmax)
	{
                y=y-Lmax;
	}
        else if (y<0)
	{
                y=y+Lmax;
	}
}

//Advective flow, with turbulence described by the Pierrehumbert map (1991)
void basic_particle::pierrehumbert_flow(double U_tot, double k, double phi, double theta, double Lmax)
{
        x=x+U_tot*cos(k*y+phi);
        y=y+U_tot*cos(k*x+theta);
        check_boundaries(Lmax);
}
