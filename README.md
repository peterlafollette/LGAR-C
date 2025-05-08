# Lumped Arid/Semi-arid Model (LASAM) for infiltration and surface runoff
The LASAM simulates infiltration and runoff based on Layered Green & Ampt with redistribution (LGAR) model. LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi-arid climates. LGAR closely mimics precipitation partitioning results simulated by the famous Richards/Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi-arid areas. LGAR in Python (no longer supported and lacking many features) is available [here](https://github.com/NOAA-OWP/LGAR-Py).

LASAM is theoretically a skillful catchment scale hydrolgic model in the event that precipitation partitioning into infiltration and runoff is the most important process for streamflow generation in a given catchment. Nonetheless, LASAM also includes a nonlinear reservoir so that some water stored in the catchment can contribute directly to streamflow (currently, we assume soil water does not directly contribute to streamflow).

A multilayer version of the Talbot-Ogden (TO) model for fluxes between groundwater and the vadose zone has been developed and coupled to LGAR, which extends the utility of the vadose zone model to humid regions with a shallow water table that does impact soil moisture dynamics. This model is called LGARTO and can be enabled with a flag in the config file. The capacity for this model to function at the catchment scale is under development (as of 7 May 2025).

Further, as part of a CIROH project, the simulation of preferential flow as dual permeability, or as simple bypass through the vadose zone, is being developed. As of this writing (7 May 2025), preferential flow as dual permeability is partially developed, and preferential flow as simple bypass is ready for use at the catchment scale. 


**Published papers:** For details about the model please see our manuscript on LGAR ([weblink](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR033742)).

## Build and Run Instructions
Detailed instructions on how to build and run LASAM can be found here [INSTALL](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/INSTALL.md).
- Test examples highlights
  - simulations with synthetic forcing data and unittest (see [build/run](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/tests/README.md)). 
  - simulations with real forcing data (see [build/run](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/INSTALL.md#standalone-mode-example))
  - LASAM coupling to Soil Freeze Thaw (SFT) model (see [instructions](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/INSTALL.md#lasam-coupling-to-soil-freeze-thaw-sft-model))

## Model Configuration File
A detailed description of the parameters for model configuration is provided [here](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/configs/README.md).

## Calibratable parameters
A detailed description of calibratable parameters is provided [here](https://github.com/peterlafollette/LGAR-C/tree/preferential_flow/data).

## Nextgen Realization Files
Realization files for running LASAM (coupled/uncoupled modes) in the nextgen framework are provided [here](https://github.com/peterlafollette/LGAR-C/tree/preferential_flow/realizations).
  
## Getting help
For questions, please contact Ahmad (ahmad.jan(at)noaa.gov) and/or Peter (peter.lafollette(at)noaa.gov), the two main developers/maintainers of the repository.

## Known issues or raise an issue
LASAM is a newly developed model and we are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/peterlafollette/LGAR-C/blob/preferential_flow/CONTRIBUTING.md)) or simply fork the repository and submit a pull request.
