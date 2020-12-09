# OneDimDVR
Numerically solve the time dependent schrodinger equation by discrete variable representation method

Analyzations of the solved wave function are available

This is a specialized version for 1 dimensional systems

## Featured utilities
`OneDimDVR-bound.exe`:
* Specialized for bound systems
* Outputs wave function

`OneDimDVR-scatter.exe`:
* Specialized for scatter systems
* Outputs wave function or transmission rate

`OneDimDVR-Chebyshev.exe`:
* Specialized for scatter systems
* Outputs Chebyshev order domain wave function

Analyzers:
* `Che2wfn.exe` converts Chebyshev order domain wave function to usual time domain
* `Che2tran.exe` calculates transmission and reflection from Chebyshev order domain wave function
* `density.exe` converts wave function to density
* `transmission.exe` calculates transmission and reflection from wave function
* `animate-density.py` animates density propagation
* `Wigner.exe` converts wave function to Wigner distribution
* `animate-Wigner.py` animates Wigner distribution propagation
* `expectation.exe` calculates expectations from Wigner distribution
* `visualize-expectation.py` visualizes expectations and makes tables for future advanced plots

## Installation
To install OneDimDVR wave function solvers, grab the `makefile` in `bound/`, `scatter/` or `Chebyshev` and see section `Link to user routines`

To install analyzers:
1. `make`, `make install`
2. All analyzers are then available in `OneDimDVR/`

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
> 3. D. E. Manolopoulos 2004 J. Chem. Phys.
> 4. V. A. Mandelshtam 1995 J. Chem. Phys.
> 5. H. Guo 2006 Phys. Rev. A