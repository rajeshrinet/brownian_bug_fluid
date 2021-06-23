//13/01/2021 CP: try and reproduce Young et al. (2001) code for the bug model

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
//gsl_rng *rgslbis2 = gsl_rng_alloc(gsl_rng_taus);

//Define constant
extern const double pi=3.14159265;
extern const double Delta=pow(10,-7); //diffusion
extern const double Lmax=pow(10,0.5); //size of the grid
extern const double area=Lmax*Lmax; //size of the grid
extern const double k=2*pi;// We do not divide by Lmax, i.e. we do not change the wavenumber even though Lmax changes 
extern const int size_pop=200000; //initial number of particles
extern const int tmax=1000; //length of the simulation
extern const double proba_death=0.5; //Death and birth probability
extern const double proba_repro=0.5;

using namespace std;

//This function computes the values of the pair density, either deriving with regards to dx, the increment in radius, or with regards to dp, the increment of power to compute the radius, as we wanted to compare the results. PairDens values are put in the table pcf[], with the derivation wrt dx in [0] and wrt to dp in [1]
double PairDens(double p, double dp, std::vector<basic_particle> Part_table,double pcf[])
{
double d2,dt2,dtt2;
int ki,kj;
basic_particle temp, current;
int iter=0;
int p1=0,p2=0;
double xi=pow(10, p);
double dxi=pow(10, p+dp)-pow(10,p);

    for(p1=0;p1<Part_table.size();p1++)
{
    current=Part_table.at(p1);
    // double loop (or loop on all pairs)
    for (p2=0;p2<Part_table.size();p2++)
{
        if(p1!=p2)
{               
                temp=Part_table.at(p2);
                dtt2=area;
                for (ki=-1;ki<=1;ki++)
            {
                    for (kj=-1;kj<=1;kj++)
                    {   
                        dt2 = pow((temp.get_x() - current.get_x()+ki*Lmax),2) + pow((temp.get_y() - current.get_y()+kj*Lmax),2);
                        d2 = std::min(dtt2,dt2);
                        dtt2 = d2;
                      }
             }

                // squared distance computed
                // if distance between xi and xi+dxi increment PairDens
                if((d2<pow(xi+dxi,2))&&(d2>pow(xi,2)))
                {
                        iter++;
                }
              }
         
        }
}
	pcf[0]=iter/(pi*(pow(xi+dxi,2.0) - pow(xi,2.0))*area);
	pcf[1]=iter/(2*pi*log(10)*pow(10,2*p))*(1/dp)*(1/area);
}


void distrib_distance(std::vector<basic_particle> Part_table, int repart[]) //This function outputs all distances between points. Should be used ONLY for small populations, and for debugging reasons.
{
double d2,dt2,dtt2;
int ki,kj;
basic_particle temp, current;
int iter=0;
int p1=0,p2=0;
int pow_dist,id_pow;
std::ofstream f3;


f3.open("distance_table.txt");

    for(p1=0;p1<Part_table.size();p1++)
{
    current=Part_table.at(p1);
    // double loop (or loop on all pairs)
    for (p2=0;p2<Part_table.size();p2++)
{
        if(p1!=p2)
{
                temp=Part_table.at(p2);
                dtt2=area;
                for (ki=-1;ki<=1;ki++)
            {
                    for (kj=-1;kj<=1;kj++)
                    {
                        dt2 = pow((temp.get_x() - current.get_x()+ki*Lmax),2) + pow((temp.get_y() - current.get_y()+kj*Lmax),2);
                        d2 = std::min(dtt2,dt2);
                        dtt2 = d2;
                      }
             }

                 if(p2<p1){
                 	pow_dist=int(std::max(int(-10),int(round(log10(pow(d2,0.5))))));
                        id_pow=-1*pow_dist;
                        repart[id_pow]=repart[id_pow]+1;
                        f3<<p1<<";"<<p2<<";"<<pow(d2,0.5)<<std::endl;
                        }

}
}
}
f3.close();
}

//Births and deaths in a community. The current vector of particles is updated with the new particles (at the end of the table), then we remove the dead ones
void branching_process(std::vector<basic_particle> &part_1,double proba_repro, double proba_death)
{
	double proba;	
	std::vector <int> id_to_remove;
	int size=part_1.size();
	basic_particle tmp_part;

	for(int j=0;j<size;j++)
	{
		proba=gsl_rng_uniform(rgslbis2); //Computes a number between 0 and 1, then compare it with the different probabilites of death, birth, or just staying alive
		if (proba<proba_repro) //new particle at the same place, and the parent is still there
		{
			tmp_part=basic_particle(part_1.at(j).get_x(),part_1.at(j).get_y(),part_1.at(j).get_yfirst(),part_1.at(j).get_firstparent());
			part_1.push_back(tmp_part);

               } else if( (proba_repro<=proba) and (proba<(proba_repro+proba_death)) )
		{
			id_to_remove.push_back(j); //Here, we just do a list of the indices we will need to remove at the end. We cannot remove them on the fly because it would modify the sequence of rows
		}
	}
	for(int j=id_to_remove.size()-1;j>=0;j--) //We start at the end of the table. Indeed, if we were to start by the beginning, removing row number 2 would shift all following rows (element 4 from the old table would be element 3, for example)
	{
			part_1.erase(part_1.begin()+id_to_remove[j]);
	}
}
 
int main()
{
	int i,j,t,tmp_t,nb_div;
	double a_x,a_y,phi,theta,a_n,xi,dxi,pow_min,pow_max,dpow,pow_i,C;
	std::vector<basic_particle> Part_table,Part_table_tmp;
	//std::vector<double> Utot_list{ 0.0, 0.1, 0.5,2.5 }; //The structure of the code enables lauching one simulation for all Utot. For speed purposes, though, it is preferable to launch one simulation per Utot
	std::vector<double> Utot_list{2.5};
	std::ofstream f0,f1;
	double pcf_table[2];

	pow_min=-1+log10(Delta);pow_max=5.5+log10(Delta); //These are the limits in Fig. 3 of Young et al. 2001
	dpow=0.25;
	int repart[11];
	for(i=0;i<11;i++){
		repart[i]=0;
	}

	//Open the file in which we will have the x, y, parent of each particle
	f0.open("nb_individuals_dpow0p25_area10_tmax1000_N200000_U2p5.txt");
	f0<<"Utot;Nb_ind;area"<<std::endl;
	f1.open("pcf_dpow0p25_area10_tmax1000_N200000_U2p5.txt");
	f1<<"r;Utot;pcf_dx;pcf_dp"<<std::endl;

	for (double Utot : Utot_list) 
	{
	//Initialize
	for(i=0; i < size_pop; i++)
	{
		a_x=gsl_rng_uniform(rgslbis2)*Lmax;
		a_y=gsl_rng_uniform(rgslbis2)*Lmax;
		Part_table.push_back(basic_particle(a_x,a_y,a_y,i));
	}

	//Run the simulation
	for(t=0;t<=tmax;t++)
	{
	std::cout<<"TIME="<<t<<std::endl;
	//Compute the phase in x and y for the turbulent flow from Pierrehumbert. These phases are common to each particle as they correspond to a unique flow
	phi=gsl_rng_uniform(rgslbis2)*2*pi;
	theta=gsl_rng_uniform(rgslbis2)*2*pi;

	//First step: birth and death	
	branching_process(Part_table,proba_repro,proba_death);

	//For each particle, diffusion, then advection
	for(j=0; j<Part_table.size(); j++)
	{
		Part_table[j].diffusion(Delta, Lmax);
		Part_table[j].pierrehumbert_flow(Utot, k, phi,theta, Lmax);

	} //end j, i.e. the particle mvt
	} //end t, i.e the whole simulation

	//End of the simulation
	pow_i=pow_min;
	C=Part_table.size()/area;
	f0<<Utot<<";"<<Part_table.size()<<";"<<area<<std::endl;
	while (pow_i<=pow_max)
	{
	
		xi=pow(10,pow_i);
		dxi=pow(10,pow_i+dpow)-pow(10,pow_i);
		f1<<xi<<";";
		f1<<Utot<<";";
		PairDens(pow_i,dpow,Part_table,pcf_table);
		f1<<pcf_table[0]/pow(C,2)<<";"<<pcf_table[1]/pow(C,2)<<std::endl;
		pow_i=pow_i+dpow;
	}

	//The part below is just here to output the repartition of distances, which should NOT be done for large populations	
	//distrib_distance(Part_table,repart);
	//for(i=0;i<11;i++){
	//	f0<<Utot<<";"<<i<<";"<<repart[i]<<std::endl;
	//	repart[i]=0; //after writing the right number in the file, we reset to 0 before a new simulation with a new Utot
	//}
	
	 Part_table= std::vector<basic_particle>(); //Deallocate, reinitialize

	} //End loop on Utot

	
	f1.close();
	f0.close();
	return 0;
}

int main();
