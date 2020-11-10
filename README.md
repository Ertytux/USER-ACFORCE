

# USER-ACFORCE
This package implements the action of standing acoustic waves by a perturbation field to be included into LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories, Steve Plimpton, sjplimp@sandia.gov

This package was created by Edisel Navas Conyedo at the University of Informatics Sciences on Cuba enavas@uci.cu, contact him directly if you have questions.

# Usage example

```
#                         xi   density  soundV     kt      kx    ky    kz  eta    beta    psi
fix  3 all  addacforce  ${xi2} ${densi} ${sounde}  ${kts} ${kx} ${ky} 0.0  ${eta} ${beta} ${psi}
```
* $$ xi^2 $$ Acoustic force strength 
* $$ densi $$ particles density 
* kx ky kz wavevector components with modulus kt
* $$ eta $$ solvent viscosity
* $$ beta $$ viscosity ratio
* $$ psi $$ pressure second order factor $$ P= densi KbT + psi densi^2 $$

