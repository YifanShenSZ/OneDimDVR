# OneDimDVR
Numerically solve the time dependent schrodinger equation by discrete variable representation method

This is a specialized version for 1 dimensional problems

## Link to user routines
OneDimDVR requires a several user routines to run:
1. The diabatic surfaces
2. The initial wave function

These routines should be provided in modules `libHd` and `libwfn` in `libHd.f90` and `libwfn.f90` with a standard interface:
* `subroutine initialize_libHd()`
* `subroutine compute_Hd(q, H)`
* `subroutine initialize_libwfn()`
* `subroutine init_wfn(q, psy)`

Linking by source (`libHd.f90` and `libwfn.f90`), static library (`libHd.a` and `libwfn.a`), dynamic library (`libHd.so` and `libwfn.so`) are all possible through setting the `link_type` variable in `makefile` to `source`, `static`, `dynamic`

The user files mentioned above should be in a same directory, whose path should be passed to `make` through the `usr_dir` variable in `makefile`

## Dependency
* My Fortran-Library, as written in `FL_dir` variable in `makefile`

## Reference
> 1. W. H. Miller 1992 J. Chem. Phys.
> 2. D. E. Manolopoulos 2002 J. Chem. Phys.