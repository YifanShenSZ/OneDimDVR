###############################################
#                                             #
#        Makefile for OneDimDVR-Wigner        #
#                                             #
###############################################

flag = -O3

OneDimDVR-Wigner.exe: main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@

.PHONY: clean
clean:
	rm -f *.mod *.o