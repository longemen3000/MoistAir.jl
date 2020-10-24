

"Molecular weight of air in kg/mol"
const Ma = 0.0289635

"Molecular weight of water in kg/mol"
const Mv = 0.01801528

"Universal gas constant kg m^2 /s^2 /K /mol"
const R = 8.314459848

"Melting point of water"
const T0 = 273.15

abstract type AbstractMoistAirModel <: ThermoState.Types.ThermoModel end
"""
Type used to model dry air

This uses the formulations presente in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5]. 
"""
struct ASHRAEDryAir <: AbstractMoistAirModel end

const DryAirModel = ASHRAEDryAir
molecular_weight(model::ASHRAEDryAir) = Ma
"""
Type used to model saturated water vapor

This uses the formulations present in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5].

"""
struct ASHRAEWaterVapor <:  AbstractMoistAirModel end

const WaterVaporModel = ASHRAEWaterVapor
molecular_weight(model::ASHRAEWaterVapor) = Mv

"""
Type used to model moist air

This uses the formulations presente in the ASHRAE handbook Psychrometrics: Theory and Practice, reference [5]. 

The temperature range is `173.15 < T < 474.15 K` with pressures `P < 5 MPa`.

It should be noted that the specific enthalpy, entropy and volume are calculated with respect to the mass of dry air. Therefore, 1 J/kg is actually 1 J per kg of dry air.
"""
struct ASHRAEMoistAir <:  AbstractMoistAirModel end
const MoistAirModel = ASHRAEMoistAir
molecular_weight(model::ASHRAEMoistAir) = SVector(1000*Ma,1000*Mv)

# package code goes here
include("utilities.jl")
include("docs.jl")
include("hyland83.jl")
include("hyland83a.jl")
include("wetbulb.jl")
include("moistair_impl.jl")


