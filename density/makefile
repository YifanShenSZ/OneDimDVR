################################################
#                                              #
#        Makefile for OneDimDVR-density        #
#                                              #
################################################

flag = -O3

OneDimDVR-density.exe: main.f90
	ifort -qopenmp -parallel -ipo $(flag) $^ -o $@

.PHONY: clean
clean:
	rm -f *.mod *.o