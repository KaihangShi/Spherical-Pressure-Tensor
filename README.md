# Spherical-Pressure-Tensor
FORTRAN code for calculating local pressure tensor in spherical coordinates. Examples here are for hard sphere fluids and TIP4P/Ice water.<br/>

![Normal and tangential pressure profile in a hard-sphere critical nucleus](https://github.com/KaihangShi/Spherical-Pressure-Tensor/blob/main/hard_sphere_fluid.jpg)


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

### Configure calculation parameters
All input parameters for the pressure tensor calculations are set at the beginning of the source  code. GROMACS trajectory file (.xtc) and coordinate file storing center-of-mass for spherical nucleus (CM_True.dat) are needed. An example of CM_true.dat is provided in the parent directory.


## Reference
[1]  P. Montero de Hijes, **K. Shi**, E. G. Noya, E. E. Santiso, K. E. Gubbins, E. Sanz and C. Vega, \"The Young–Laplace equation for a solid–liquid interface\". *Journal of Chemical Physics*, 153 (2020) 191102. [[link]](https://aip.scitation.org/doi/10.1063/5.0032602)[[PDF]](http://kaihangshi.github.io/assets/docs/paper/Hijes_jcp_2020.pdf)
