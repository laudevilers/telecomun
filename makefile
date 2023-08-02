# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -Wall -g
 
# ****************************************************
# Targets needed to bring the executable up to date
 
main: main.o norme.o vecteur.o coefficient_reflexion.o gamma_m.o t_m.o mur.o comp_directe.o
	$(CC) $(CFLAGS) -o main main.o norme.o vecteur.o coefficient_reflexion.o gamma_m.o t_m.o mur.o comp_directe.o
 
# The main.o target can be written more simply
 
main.o: main.cpp norme.cpp vecteur.hpp coefficient_reflexion.hpp gamma_m.hpp t_m.hpp mur.hpp comp_directe.hpp
	$(CC) $(CFLAGS) -c main.cpp
	
norme.o: norme.hpp
 
vecteur.o: vecteur.hpp

coefficient_reflexion.o: coefficient_reflexion.hpp

gamma_m.o : gamma_m.hpp

t_m.o : t_m.hpp

mur.o :  mur.hpp

comp_directe.o: comp_directe.hpp
