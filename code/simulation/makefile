main_Fig3.out : main_Fig3.o basic_particle.o 
	g++ -O3 -std=c++11 -o main_Fig3.out basic_particle.o main_Fig3.o -lgsl -lgslcblas -lm
main_Fig3.o : main_Fig3.cpp basic_particle.h
	g++ -O3 -std=c++11 -o main_Fig3.o -c main_Fig3.cpp
basic_particle.o : basic_particle.h
	g++ -c -Wall -std=c++11 -o basic_particle.o basic_particle.cpp -lgsl -lgslcblas -lm

