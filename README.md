# Spherical-Pressure-Tensor
FORTRAN code for calculating local pressure tensor in spherical coordinates. <br/>

![Normal and tangential pressure profile in a hard-sphere critical nucleus](https://github.com/KaihangShi/Spherical-Pressure-Tensor/blob/main/hard_sphere_fluid.jpg)


## Files
- Foler `/TIP4P_Ice/` contains source code to calculate spherical and planar local pressure tensor profile using TIP4P/Ice model, as described in Ref. 1. <br/>
- Folder `/Hard_sphere/` contains source code used to calculate spherical pressure tensor using hard sphere model, as described in Ref. 2. 

## Requirement
1). gfortran compiler <br/>
2) [libgmxfort](https://github.com/KaihangShi/libgmxfort) for reading GROMACS compressed trajectory files (.xtc)<br/>

## Usage
### Compilation
Once installed libgmxfort, FORTRAN code for pressure tensor calculation can be compiled using
```bash
gfortran $file -I/usr/local/include `pkg-config --libs libgmxfort`

```
where ```/usr/local/include``` is path to the installed module file. Replace ```$file``` with the corresponding file name for FORTRAN code.

### Computation parameters
All input parameters for the pressure tensor calculations are set at the beginning of the FORTRAN source code. GROMACS trajectory file (.xtc) and center-of-mass coordinate file for spherical nucleus (CM_True.dat) are needed. "CM_True.dat" file is used to translate the system origin to the center of the nucleus. An example of "CM_true.dat" is provided in the `/Hard_sphere/` directory.

## Reference
[1] P. Montero de Hijes=, **K. Shi=**, C. Vega, and C. Dellago, \"Comparing the Mechanical and Thermodynamic Definitions of Pressure in Ice Nucleation\". *to be submitted*. <br/>
[2] P. Montero de Hijes, **K. Shi**, E. G. Noya, E. E. Santiso, K. E. Gubbins, E. Sanz and C. Vega, \"The Young–Laplace equation for a solid–liquid interface\". *Journal of Chemical Physics*, 153 (2020) 191102. [[link]](https://aip.scitation.org/doi/10.1063/5.0032602)[[PDF]](http://kaihangshi.github.io/assets/docs/paper/Hijes_jcp_2020.pdf)
