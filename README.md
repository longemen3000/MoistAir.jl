# MoistAir - Thermodynamic properties of moist air

[![Build Status](https://github.com/longemen3000/MoistAir.jl/workflows/CI/badge.svg)](https://github.com/longemen3000/MoistAir.jl/actions)
[![Build Status](https://travis-ci.com/longemen3000/MoistAir.jl.svg?branch=master)](https://travis-ci.com/longemen3000/MoistAir.jl)
[![Codecov](https://codecov.io/gh/longemen3000/MoistAir.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/longemen3000/MoistAir.jl)

This package provides a model, `MoistAirModel` to compute thermodynamic properties of moist air, using the [ThermoState](https://github.com/longemen3000/ThermoState.jl) interface. The model uses real gas correlations as recommended by ASHRAE (see reference [5]). 

This package also provides the underlying models used for moist air, namely, dry air (`ASHRAEDryAir`) and water vapor (`ASHRAEWaterVapor`). those models have the same interface as the moist air model.

this is a fork of [Psychro.jl](https://github.com/pjabardo/Psychro.jl/) working with julia 1.5+.

## instalation

```julia-repl
julia> ]
(v1.5) pkg> add https://github.com/longemen3000/MoistAir.jl
```
## User interface - Thermodynamic properties of moist air, dry air and saturated water vapor.

The models defined in this package accept the following thermodynamic states:
- `ASHRAEDryAir`: pressure,temperature
- `ASHRAEWaterVapor`: pressure,temperature
- `MoistAirModel` : 
    - pressure, temperature, amounts (any of mol_fraction,mass_fraction,mol_number,mass_number)
    - pressure,temperature, humidity spec


The models defined in this package support the following `ThermoState` functions:

- `mol_enthalpy`, `mass_enthalpy`, `total_enthalpy` 
- `mol_entropy`, `mass_entropy`, `total_entropy` 
- `mol_volume`, `mass_volume`, `total_volume`
- `mol_enthalpy`, `mass_enthalpy`
- `compressibility_factor`

The package also provides its own functions, with ThermoState syntax, specific for moist air:

- Dry properties: properties per 1 kg of dry air:
    - `dry_volume`,`dry_enthalpy`,`dry_entropy`

- `rel_hum` : Relative humidity, defined as the fraction between the amount of water in the air and the amount of water of saturated air at the same temperature and pressure (between 0 and 1)

- `hum_massfrac` : kg of water / kg of moist air

- `hum_molfrac` : mol of water / mol of moist air

- `hum_wetbulb` : wet-bulb temperature or adiabatic saturation temperature is the temperature a volume of air would have if cooled adiabatically to saturation by evaporation of water into it, always lower than the specified temperature (also known as dry-bulb temperature)

- `hum_ratio` : Humidity ratio (kg of vapor / kg of dry air)

- `hum_dewpoint`  temperature to which moist air must be cooled to become saturated,always lower than the specified temperature
## Humidity Spec

As mentioned before, the MoistAirModel accepts thermodynamic states of the form pressure-temperature-humidity spec. you can use the same function names as keywords in the `ThermoState.state` function.

## Usage

```julia
using ThermoState,Unitful,MoistAir
air = MoistAirModel() #Moist air model
λ = VariableSpec() #to create a variable model
st = state(hum_wetbulb = 25.0u"°C",t=λ,p=1u"atm")
dryh = dry_enthalpy(air,st(32u"°C")) #custom units
x = mol_fraction(air,st(300.15)) #molar fraction of the mixture air-water
xw = hum_molfrac(air,st(300.15)) #molar fraction of water
dew = hum_dewpoint(air,st(300.15)) #dew point of water
last(x) == xw #true
```

## References

 * [1] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of the saturated phases of H2O from 173.15 K to 473.15 K", ASHRAE Transactions, 1983.
 * [2] Wexler, A. and Hyland, R. W., "Formulations for the thermodynamic properties of dry air from 173.15 K to 473.15 K, and of saturated moist air from 173.15 K to 372.15 K at pressures to 5 MPa
 * [3] Himmelblaum D. M., "Solubilities of inert gases in water, 0oC to near the critical point of water", Journal of Chemical and Engineering Data, Vol. 5, No. 1, January 1960.
 * [4] Kell, George S., "Density, thermal expansivity, and compressibility of liquid water from 0oC to 150oC: correlations and tables for atmospheric pressure and saturation reviewed and expressed on 1968 temperature scale", Journal of Chemical and Engineering Data, Vol. 20, No. 1, 1975.
 * [5] ASHRAE, "Psychrometrics: Theory and Practice", ASHRAE, 1996.
