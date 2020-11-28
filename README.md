# OneDimDVR
Numerically solve the time dependent schrodinger equation by discrete variable representation method

This is a specialized version for 1 dimensional systems

## Featured utilities
bound: `OneDimDVR-boud.exe`
* Specialized for bound systems
* Outputs wave function

scatter: `OneDimDVR-scatter.exe`
* Specialized for scatter systems
* Outputs wave function or transmission & reflection

density: `density.exe` and `animate.py`
* Convert wave function to density
* Animates density propagation

Wigner: `Wigner.exe` and `animate.py`
* Convert wave function to Wigner distribution
* Animates Wigner distribution propagation

## Link to user routines
OneDimDVR requires a several user routines to run:
1. The number of electronic states
2. The diabatic surfaces
3. The initial wave function

These routines should be provided in modules `libHd` and `libwfn` with a standard interface:
1. module libHd
* `integer, parameter::NStates =`
* `subroutine initialize_libHd()`
* `subroutine compute_Hd(q, H)`
2. module libwfn
* `subroutine initialize_libwfn()`
* `subroutine init_wfn(q, psy)`

Linking by source (`libHd.f90` and `libwfn.f90`), static library (`libHd.a` and `libwfn.a`), dynamic library (`libHd.so` and `libwfn.so`) are all possible through setting the `link_type` variable in `makefile` to `source`, `static`, `dynamic`

The user files mentioned above should be in a same directory, whose path should be passed to `make` through the `usr_dir` variable in `makefile`

## Dependency
* My Fortran-Library, as written in `FL_dir` variable in `makefile`

## Reference
> 1. W. H. Miller 1992 J. Chem. Phys.
> 2. D. E. Manolopoulos 2002 J. Chem. Phys.