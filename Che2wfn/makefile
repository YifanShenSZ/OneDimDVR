######################################
#                                    #
#        Makefile for Che2wfn        #
#                                    #
######################################

flag = -O3

Che2wfn.exe: main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@
