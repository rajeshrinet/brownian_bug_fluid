main_Fig1.out : main_Fig1.o basic_particle.o 
	g++ -std=c++11 -o main_Fig1.out basic_particle.o main_Fig1.o -lgsl -lgslcblas -lm
main_Fig1.o : main_Fig1.cpp basic_particle.h
	g++ -std=c++11 -o main_Fig1.o -c main_Fig1.cpp
basic_particle.o : basic_particle.h
	g++ -c -Wall -std=c++11 -o basic_particle.o basic_particle.cpp -lgsl -lgslcblas -lm

