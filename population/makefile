#########################################
#                                       #
#        Makefile for population        #
#                                       #
#########################################

flag = -O3

population.exe: main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@
