cd ../simulation 
echo simulating using code in the simulation folder 
echo -------------------------------------
echo 
g++ -c -Wall -std=c++11 -o basic_particle.o basic_particle.cpp -lgsl -lgslcblas -lm
g++ -c -Wall -std=c++11 -o main_Fig1.o main_Fig1.cpp -lgsl -lgslcblas -lm
g++ -O3 -std=c++11 -o main_Fig1.out basic_particle.o main_Fig1.o -lgsl -lgslcblas -lm 
./main_Fig1.out 

cd ../figure/
echo 
echo -------------------------------------
echo plotting using code in the simulation folder 
echo -------------------------------------
echo 
Rscript visualisation_Fig1.r
