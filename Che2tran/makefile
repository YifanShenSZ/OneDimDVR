#######################################
#                                     #
#        Makefile for Che2tran        #
#                                     #
#######################################

flag = -O3

Che2tran.exe: main.f90
	ifort -qopenmp -parallel -ipo -mkl $(flag) $^ -o $@
