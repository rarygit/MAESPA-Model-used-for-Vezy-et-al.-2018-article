# MAESPA 0.1.1

* Added a `NEWS.md` file to track changes.
* Add semantic versioning as of [semver specifications](https://semver.org/) of the form `MAJOR.MINOR.PATCH`
* The [Makefile](Makefile) is now compatible with Windows OS.
* Fix error in parameter reading for the Medlyn stomatal conductance model: `G0` and `G1` were expected for `CO2`, but the model said it was expected for `H2O`, and it did not perform any conversion as for the other models.
* `ET` and `EV` can be negative when there is dew. If `ET` is negative, we now transfer its values to `EV` to increase CANOPY_STORE for the next time-step, and set it equal to `0.0` (see commit f4fda2faea2792c5bfd141b41bfb1a22da497302).
* Re-activate the initialization values in the `watbal.dat` (can be informative for the user, but mostly because `Maeswrap::readwatbal()` handles only the file when this line is present.  
* Add a [logo](logo/MAESPA.png)

# MAESPA 0.1.0

First release of the code based on the [`Vezy_et_al_2018`](https://bitbucket.org/remkoduursma/maespa/branch/Vezy_et_al_2018) branch from the historic repository from R. Duursma and B. Medlyn that already implemented the following changes:

- voxel-scale leaf evaporation (EV)
- Computation of canopy air temperature and vapor pressure.
- Iterations on leaf temperature until convergence
- Aerodynamic conductance computation at plot scale (used for soil + canopy air temperature and vapor pressure)
- Allowing two wind extinction profile if ZHT (measurement height) is below canopy height (mainly used in agroforestry with measurements below the shade trees)
- Bug fixes on slope computation

The current master branch of MAESPA in this repository as of version `0.1.0` is then based on the [`Vezy_et_al_2018`](https://bitbucket.org/remkoduursma/maespa/branch/Vezy_et_al_2018) branch of the original Bitbucket project, implemented with the following small modifications by Céline Blitz Frayret (modifications in the source code are tagged by 'CBF'):

- Add of run timing screen prints
- Modifications of the static allocations of the following variables in files maindeclarations.f90 and watbal.f90 : wettingbot, wettingtop, plantwater
- Modification of the value of maxdate parameter in maestcom.f90
- Inactivation of a section of watbal.f90 to prevent from dummy initial values in the output *_watbal.dat file
