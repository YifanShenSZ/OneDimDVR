# OneDimDVR
Numerically solve the time dependent schrodinger equation by discrete variable representation method
This is a special version for 1 dimensional problem

Supported job types:
1. NewTrajectory
* Propagate wavefunction based on user specified initial condition, and save the trajectory
2. TR-p0
* Compute transmission and reflection on each state respect to different initial momentum
3. SMD
* Calculate SMD quantities based on an old trajectory
4. pRepresentation
* Convert an old trajectory (in coordinate representation) to momemtum representation
5. WignerDistribution
* Convert an old trajectory to phase space

Analyze.py provides visualization for each type of job. Just run it directly after OneDimDVR.exe finishes

Modify Basic.f90 for different initial wave function and potential

Dependency:
* My Fortran-Library, as written in MyLib of makefile
* My Python-Library, as written in Import library section of Analyze.py

Reference:
> 1. W. H. Miller 1992 J. Chem. Phys.
> 2. D. E. Manolopoulos 2002 J. Chem. Phys.