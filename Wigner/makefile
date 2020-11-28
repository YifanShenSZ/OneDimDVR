#####################################
#                                   #
#        Makefile for Wigner        #
#                                   #
#####################################

flag = -O3

Wigner.exe: main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@

.PHONY: clean
clean:
	rm -f *.mod *.o