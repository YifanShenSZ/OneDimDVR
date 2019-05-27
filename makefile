########################################
#                                      #
#        Makefile for OneDimDVR        #
#                                      #
########################################

compiler = ifort
MyLibDir = /home-4/yshen57@jhu.edu/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90
src = Basic.f90 DVR.f90 Main.f90
exe = OneDimDVR.exe
flag = -u -mkl -fast -march=core-avx2

$(exe): $(MyLib) $(src)
	$(compiler) $(flag) $^ -o $(exe)

clean:
	rm *.mod