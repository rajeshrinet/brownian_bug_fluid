cd ../simulation 
echo simulating using code in the simulation folder 
echo -------------------------------------
echo 
g++ -c -Wall -std=c++11 -o basic_particle.o basic_particle.cpp -lgsl -lgslcblas -lm
g++ -c -Wall -std=c++11 -o main_Fig2a.o main_Fig2a.cpp -lgsl -lgslcblas -lm
g++ -O3 -std=c++11 -o main_Fig2a.out basic_particle.o main_Fig2a.o -lgsl -lgslcblas -lm 
./main_Fig2a.out 

cd ../figure/
echo 
echo -------------------------------------
echo plotting using code in the simulation folder 
echo -------------------------------------
echo 
Rscript visualisation_Fig2.r
