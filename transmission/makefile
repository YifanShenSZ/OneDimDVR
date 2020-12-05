###########################################
#                                         #
#        Makefile for transmission        #
#                                         #
###########################################

flag = -O3

transmission.exe: main.f90
	ifort -qopenmp -parallel -ipo -mkl $(flag) $^ -o $@
