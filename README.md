# MAESPA
A soil-plant-atmosphere model based on MAESTRA and SPA <img src="logo/MAESPA.png" alt="logo" width="300" align="right" />

---
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

Maespa is a model of forest canopy radiation absorption, photosynthesis and water balance. The model has a long history, going back to the work of John Norman and Paul Jarvis in the 1970's and 80's. Ying-Ping Wang improved and tested the model for his PhD thesis (then called MAESTRO) (Wang & Jarvis 1990). In 1997 Belinda Medlyn obtained a version of the code from Ying-Ping Wang and revised the model, with the objectives of:
1. modularising the code to make the program easier to understand and modify; and
1. incorporating standard formulations of leaf gas exchange models (see Medlyn 2004 for a full history).

In 2008, Remko Duursma and Belinda Medlyn started expanding Maestra by including water balance and plant hydraulics routines, largely based on the SPA model by Mat Williams (Duursma and Medlyn 2012). The two models were independent for some time, but since 2013 they use the same code base. It is straightforward to run the model 'in Maestra mode', which is consistent with the legacy Maestra model. MAESPA is no longer maintained by R. Duursma and B. Medlyn, as explained on its [original GitHub website](https://maespa.github.io/).

So as to pursue developments and use of MAESPA, the administration of a particular branch ([`Vezy_et_al_2018`](https://bitbucket.org/remkoduursma/maespa/branch/Vezy_et_al_2018) developed by Rémi Vezy during its PhD at CIRAD) has been attributed to CIRAD, under common agreements, and is hosted on this repository.

Please take a look at the [news](NEWS.md) file for changes.

## R package

A companion R package is distributed on [CRAN](https://cran.r-project.org/web/packages/Maeswrap/index.html) and on [Github]([Maeswrap](https://github.com/RemkoDuursma/Maeswrap/tree/master/R)).

Download the last CRAN version using:

```r
install.packages("Maeswrap")
```

Or the latest development version using:

```r
devtools::install_github("RemkoDuursma/Maeswrap")
# Or if you don't have devtools:
remotes::install_github("RemkoDuursma/Maeswrap")
```

## Compilation

Compiles with:
* gfortran
* Intel Visual Fortran Compiler (version >10).

A Makefile is provided to compile MAESPA on any platform using `gfortran`. To use it, open a terminal in the root folder, and type one of these commands:

 * `make all`: build maespa in the build directory (must already exist);
 * `make clean`: remove the files in the build directory.

Windows users can also build maespa using `make` and `gfortran` by following the steps described in the file [`windows_gfortran_info.md`](windows_gfortran_info.md).

## Further information

Please note that this version of MAESPA is suitable for a "distributed" run (further sequential runs launched at the same time), so input parameters are now numbered (e.g. 1_confile.dat for the first plot, 2_confile.dat for the 2nd plot if any, etc …). Once the executable of MAESPA is activated, a number corresponding to the simulated plot must be provided.


This MAESPA project is released under the same license of the original Bitbucket project of MAESPA : the GNU General Public License Version 2.

## Contributors
MAESPA creators: [Remko Duursma](http://www.remkoduursma.com/) and [Belinda Medlyn](https://bmedlyn.wordpress.com) developed MAESPA from `MAESTRA` and `SPA`.  

Authors (alphabetical order): [Céline Blitz-Frayret](https://umr-ecosols.fr/en/eco-sols/annuaire/fiches-annuaire/177-blitz-frayret), [Mathias Christina](https://www.researchgate.net/profile/Mathias_Christina), [Guerric le Maire](https://www.researchgate.net/profile/Guerric_Le_Maire) and [Rémi Vezy](https://remi-vezy.netlify.com/).  

## References
Duursma RA & Medlyn BE (2012) MAESPA: A model to study interactions between water limitation, environmental drivers and vegetation function at tree and stand levels, with an example application to [CO2] x drought interactions. Geoscientific Model Development 5 : 919-940 [Link](http://www.geosci-model-dev.net/5/919/2012/gmd-5-919-2012.html)

Vezy, R., et al., Measuring and modelling energy partitioning in canopies of varying complexity using MAESPA model. Agricultural and Forest Meteorology, 2018. 253–254: p. 203-217. [Official link](https://www.sciencedirect.com/science/article/pii/S016819231830042X) or
[Researchgate](https://www.researchgate.net/publication/323398728_Measuring_and_modelling_energy_partitioning_in_canopies_of_varying_complexity_using_MAESPA_model) alternatively.
