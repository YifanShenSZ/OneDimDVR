##########################################
#                                        #
#        Makefile for expectation        #
#                                        #
##########################################

flag = -O3

expectation.exe: StringUtility.f90 polynomial.f90 main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@
