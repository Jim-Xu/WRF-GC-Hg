# WRF-GC-Hg

Here we did some modifications for `WRF-GC` to contain elaborate mercury simulation. This repository (`WRF-GC-Hg`) contains three parts, `chem` is a supplementary directory based on `WRF-GC v1.0` model (https://doi.org/10.5281/zenodo.3550330; Lin et al., 2019), `for input` contains modifications when set up simulations in `run` folder, `for ICBC` contains modifications for initial and boundary conditions.

For model set-up and installation issues, please see (https://wrf.geos-chem.org) for more instructions.
For initial and boundary condition issues, please see repositories from `Xu Feng` [GitHub Pages](https://github.com/fengx7).

Below is the path for files that have modified.
1. WRF-GC/WRFV3/`registry.chem`
2. WRF-GC/WRFV3/chem/`gigc_convert_state_mod.F`
3. WRF-GC/WRFV3/chem/gigc/KPP/Standard/`Standard.eqn`
4. WRF-GC/WRFV3/chem/gigc/KPP/Standard/`gckpp.kpp`
5. WRF-GC/WRFV3/chem/gigc/KPP/Standard/`gckpp_HetRates.F90`
6. WRF-GC/WRFV3/chem/gigc/GeosCore/`part_mercury_mod.F`
7. WRF-GC/WRFV3/chem/gigc/GCHP/`gigc_chunk_mod.F90`
8. WRF-GC/WRFV3/run/`input.geos`
9. WRF-GC/WRFV3/run/`HEMCO_Config.rc`
10. WRF-GC/WRFV3/run/`namelist.inpt`
11. WRF-GC/mozbc/`GEOS-Chem.inp`

If you have any question related to `WRF-GC`, please contact: Haipeng Lin `hplin@seas.harvard.edu`, Xu Feng `fengx7@pku.edu.cn`, Tzung-May Fu `fuzm@sustech.edu.cn`
If you have any question related to `WRF-GC-Hg`, please contact: Xiaotian Xu `xx24@illinois.edu`, Yanxu Zhang `zhangyx@nju.edu.cn`
