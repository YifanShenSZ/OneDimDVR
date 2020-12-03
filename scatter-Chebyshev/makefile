##########################################################
#                                                        #
#        Makefile for OneDimDVR-scatter-Chebyshev        #
#                                                        #
##########################################################

# Flags for user to tune
source_dir = ~/Software/Mine/OneDimDVR/scatter-Chebyshev/
FL_dir = ~/Library/Fortran-Library/
usr_dir = 
link_type = source
flag = -O3

# User does not have to take care of following variables
src = $(addprefix $(source_dir), basic.f90 DVR.f90 solver.f90 main.f90)
inc = $(addprefix $(FL_dir), include)
lib = $(addprefix $(FL_dir), lib/libFL.a)

OneDimDVR-scatter-Chebyshev.exe: $(src)
ifeq ($(link_type),source)
	ifort -qopenmp -parallel -mkl -static-intel -ipo -I$(inc) $(flag) $(addprefix $(usr_dir), libHd.f90 libwfn.f90) $^ $(lib) -o $@
else
ifeq ($(link_type),static)
	lib = $(lib) $(addprefix $(usr_dir), libHd.a libwfn.a)
	ifort -qopenmp -parallel -mkl -static-intel -ipo -I$(inc) $(flag) $^ $(lib) -o $@
else
	ifort -qopenmp -parallel -mkl -static-intel -ipo -I$(inc) $(flag) $^ $(lib) -lHd -lwfn -o $@
endif
endif

.PHONY: clean
clean:
	rm -f *.mod *.o