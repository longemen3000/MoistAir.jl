module MoistAir

    using Unitful
    using ThermoState
    using Roots

    using StaticArrays
    using ThermoState.Types
    using ThermoState.QuickStates

    import ThermoState: pressure,temperature,mass,moles,molar_mass
    import ThermoState: mass_volume, mol_volume, total_volume
    import ThermoState: mass_density, mol_density
    import ThermoState: mass_enthalpy, mol_enthalpy,total_enthalpy
    import ThermoState: mass_gibbs, mol_gibbs, total_gibbs
    import ThermoState: mass_helmholtz, mol_helmholtz, total_helmholtz
    import ThermoState: mass_internal_energy, mol_internal_energy, total_internal_energy
    import ThermoState: mass_entropy, mol_entropy, total_entropy
    import ThermoState: mass_fraction, mol_fraction
    import ThermoState: mass_number, mol_number

    include("types.jl")
#=
,WebBulbPT  => :hum_wetbulb
,HumidityRatioPT => :hum_ratio
,MolHumidityPT => :hum_molfrac 
,MassHumidityPT => :hum_massfrac
,RelativeHumidityPT => :rel_hum
,DewPointPT => :hum_dewpoint

=#
    export ASHRAEMoistAir, ASHRAEDryAir, ASHRAEWaterVapor
    export MoistAirModel

    #dry quantities
    export dry_volume,dry_enthalpy,dry_entropy,molecular_weight

    #specific humidity properties
    export hum_wetbulb,hum_dewpoint
    export hum_molfrac,hum_massfrac
    export rel_hum,hum_ratio

end