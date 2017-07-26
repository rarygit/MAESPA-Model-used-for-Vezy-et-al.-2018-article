A soil-plant-atmosphere model based on MAESTRA and SPA
===================================================

Version derived from the original [MAESPA model](https://maespa.github.io/).  
Changes:
* voxel-scale leaf evaporation
* Computation of canopy air temperature and vapor pressure.
* Iterations on leaf temperature until convergence
* Bug fixes

Original model maintainer: [Remko Duursma](http://www.remkoduursma.com/).  
Compiles with:
* Intel Visual Fortran Compiler (version >10).
* gfortran

NB: The Intel Visual Fortran Compiler can lead to much faster executables using optimisations.  
A Makefile is provided to compile Maes* on a Mac or linux-based OS (thanks to Martin de Kauwe and Alejandro Morales).

Here are the links to direct download this version executables for Windows [Win x84](https://gitlab.com/VEZY/MAESPA/raw/3bcfa81e53840c609f50a7a5d948257cc33a0822/Executables/Windows/maespa_x86.exe) and [Win x64 ](https://gitlab.com/VEZY/MAESPA/raw/3bcfa81e53840c609f50a7a5d948257cc33a0822/Executables/Windows/maespa_x64.exe) and Linux [CentOS](https://gitlab.com/VEZY/MAESPA/raw/3bcfa81e53840c609f50a7a5d948257cc33a0822/Executables/CentOS/maespa), [Debian](https://gitlab.com/VEZY/MAESPA/raw/3bcfa81e53840c609f50a7a5d948257cc33a0822/Executables/Debian/maespa) and [RedHat](https://gitlab.com/VEZY/MAESPA/raw/3bcfa81e53840c609f50a7a5d948257cc33a0822/Executables/RedHat/maespa).

Here the links to the official depository on which this version is based on [here](https://bitbucket.org/remkoduursma/maespa/commits/b7440289f6801fbd8841241e48ead6592e58e39e?at=montpellier) and direct download [link](https://bitbucket.org/remkoduursma/maespa/get/montpellier.zip).

To adapt the input files from the original version to this version, please consider reading the [wiki](https://gitlab.com/VEZY/MAESPA/wikis/home).
