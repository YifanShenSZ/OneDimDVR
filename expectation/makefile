##########################################
#                                        #
#        Makefile for expectation        #
#                                        #
##########################################

flag = -O3

expectation.exe: polynomial.f90 main.f90
	ifort -qopenmp -parallel -ipo -lFL $(flag) $^ -o $@
