
function dryair_check(p,t)
    if !(173.1 < t < 473.2)
        return DomainError("temperature should be between 173.15 K and 473.15 K")
    elseif !(0.0 <= p < 5e6)
        return DomainError("Pressure should be below 5×10^5 Pa")
    else
        return nothing
    end
end
function mass_volume_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)   
    return volumeair(t, p)
end


function mol_volume_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)
    return molarvolumeair(t, p)
end
    

function mass_density_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    return one(p)/mass_volume_impl(mt,model,p,t)
end

function mol_density_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    return one(p)/mol_volume_impl(mt,model,p,t)
end


function mass_enthalpy_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)    
    return enthalpyair(t, p)
end


function mol_enthalpy_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)    
    return molarenthalpyair(t, p)
end


function mass_entropy_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)    
    return entropyair(t, p)
end


function mol_entropy_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)    
    return molarentropyair(t, p)
end


function compressibility_factor_impl(mt::SinglePT,model::ASHRAEDryAir,p, t)
    dryair_check(p,t)    
    return p*molarvolumeair(t, p) / (R*t)
end


function water_vapor_check(t)
    if !(173.1 < t < 473.2)
        return DomainError("temperature should be between 173.15 K and 473.15 K")
    else
        return nothing
    end
end
function mass_volume_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return volumevapor(t)
end
    
function mol_volume_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return molarvolvapor(t)
end

function mass_density_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return one(t)/volumevapor(t)
end

function mol_density_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return one(t)/molarvolvapor(t)
end

function mass_enthalpy_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    enthalpyvapor(t)
end

function mol_enthalpy_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return enthalpyvapor(t)*Mv
end

function mol_entropy_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return molarentropyvapor(t)
end

function mass_entropy_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return entropyvapor(t)
end

function compressibility_factor_impl(mt::SingleSatT, model::ASHRAEWaterVapor,t)
    water_vapor_check(t)
    return Zvapor(t)
end

#=
,WebBulbPT  => :hum_wetbulb
,HumidityRatioPT => :hum_ratio
,MolHumidityPT => :hum_molfrac 
,MassHumidityPT => :hum_massfrac
,RelativeHumidityPT => :rel_hum
,DewPointPT => :hum_dewpoint

=#

const HumPT = Tuple{Pressure,Temperature,HumiditySpec}
const MolHumidityPT = Tuple{Pressure,Temperature,HumiditySpec{MolarHumidity}}
const MassHumidityPT = Tuple{Pressure,Temperature,HumiditySpec{MassHumidity}}
const HumidityRatioPT = Tuple{Pressure,Temperature,HumiditySpec{HumidityRatio}}
const RelativeHumidityPT = Tuple{Pressure,Temperature,HumiditySpec{RelativeHumidity}}
const DewPointPT = Tuple{Pressure,Temperature,HumiditySpec{HumidityDewPoint}}
const WebBulbPT = Tuple{Pressure,Temperature,HumiditySpec{WetBulbTemperature}}

"mixture molecular weight of moist air"
moist_mw(xv) = xv*Mv + (1-xv)*Ma
#for interaction with ThermoState
function hum_molfrac_impl(mt::MultiPT,model::ASHRAEMoistAir,p,t,x)
    return last(x)
end

function hum_molfrac_impl(mt::MolHumidityPT,model::ASHRAEMoistAir,p,t,xv)
    return xv
end

function hum_molfrac_impl(mt::HumidityRatioPT,model::ASHRAEMoistAir, p,t,hum)
    return hum  / (Mv/Ma + hum)
end

function hum_molfrac_impl(mt::MassHumidityPT,model::ASHRAEMoistAir,p,t,hum)
    return hum * Ma / (Mv + hum*(Ma - Mv))
end

function hum_molfrac_impl(mt::RelativeHumidityPT,model::ASHRAEMoistAir,p,t,hum)
    return hum * efactor(t, p) * Pws(t) / p
end

function hum_molfrac_impl(mt::DewPointPT,model::ASHRAEMoistAir,p,t,hum)
    return efactor(hum,p) * Pws(hum) / p
end

function hum_molfrac_impl(mt::WebBulbPT,model::ASHRAEMoistAir,p,t,B)
    w = calc_W_from_B(t, B, p)
    return w / (Mv/Ma + w)
end


function dry_volume_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return volumemoist(t, p, xv)
end

function mass_volume_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    mw = moist_mw(xv)
    vm = molarvolumemoist(t, p, xv)
    return vm/mw
end


function mol_volume_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return molarvolumemoist(t, p, xv)
end

function mol_density_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    vm =  molarvolumemoist(t, p, xv)
    return one(vm)/vm
end


function mass_density_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    mw = moist_mw(xv)
    vm = molarvolumemoist(t, p, xv)
    return mw/vm
end


function dry_enthalpy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return enthalpymoist(t, p, xv)
end

function mol_enthalpy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return molarenthalpymoist(t, p, xv)
end

function mass_enthalpy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    mw = moist_mw(xv)
    hn = molarenthalpymoist(t, p, xv)
    return hn/mw
end

function dry_entropy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return entropymoist(t, p, xv)
end

function mass_entropy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    mw = moist_mw(xv)
    sn = molarentropymoist(t, p, xv)
    return sn/mw
end

function mol_entropy_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return molarentropymoist(t, p, xv)
end

function compressibility_factor_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return Zmoist(t, p, xv)
end

function calcdewpoint(Tk, P, xv, EPS=1e-9, MAXITER=100)
    # Use Ideal Gas to 
    D = Tws(xv*P)
    Dnew = D
    i = 0
    err = 0.0
    for i = 1:MAXITER
        f = efactor(D, P)
        Dnew = Tws(xv*P/f)
        err = abs(D-Dnew)
        if err < EPS
            return Dnew
        end
        D = Dnew
    end
    return Dnew
end


function hum_dewpoint_impl(mt::DewPointPT, model::ASHRAEMoistAir,p,t,hum)
    return hum
end

function hum_dewpoint_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    D = calcdewpoint(t, p, xv)
    return D
end

function rel_hum_impl(mt::RelativeHumidityPT, model::ASHRAEMoistAir,p,t,hum)
    return hum
end

function rel_hum_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)   
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    xv * P / (efactor(t, p) * Pws(t))
end


function hum_wetbulb_impl(mt::WebBulbPT, model::ASHRAEMoistAir,p,t,hum)
    return hum
end
function hum_wetbulb_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)   
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    B = calcwetbulb(t, p, xv)
    return B 
end



_humrat(xv) = xv / (1-xv) * (Mv/Ma)

function hum_ratio_impl(mt::HumidityRatioPT, model::ASHRAEMoistAir,p,t,hum)
    return hum
end
function hum_ratio_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)   
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return _humrat(xv)
end


function hum_massfrac_impl(mt::MassHumidityPT, model::ASHRAEMoistAir,p,t,hum)
    return hum
end

function hum_massfrac_impl(mt::HumPT, model::ASHRAEMoistAir,p,t,hum)
    dryair_check(p,t)   
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return xv*Mv / ( (1-xv)*Ma + xv*Mv )
end

#implementations done, now with thermostate conections:

#extracts the value of the humidity spec
function humidity_spec(::FromState,st::ThermodynamicState)
    return value(get_spec(HumiditySpec,st))
end

function mol_fraction_impl(mt::HumPT,model::ASHRAEMoistAir,p,t,hum)
    xv = hum_molfrac_impl(mt,model,p,t,hum)
    return SVector(one(xv)-xv,xv)
end

function mass_fraction_impl(mt::HumPT,model::ASHRAEMoistAir,p,t,hum)
    xv = hum_massfrac_impl(mt,model,p,t,hum)
    return SVector(one(xv)-xv,xv)
end

#hum functions
for (hum_op,humunit) in zip(
    [:hum_dewpoint,:hum_massfrac,
    :hum_molfrac,:rel_hum,
    :hum_wetbulb,:hum_ratio,
    :dry_volume,:dry_enthalpy,:dry_entropy],
    [u"K",Unitful.NoUnits,
    Unitful.NoUnits,Unitful.NoUnits,
    u"K",Unitful.NoUnits,
    u"m^3/kg",u"J/kg",u"J/(kg*K)"])
    hum_op_impl = Symbol(hum_op,:_impl)
    
    @eval begin
        
        ThermoState.default_units(x::typeof($hum_op)) = $humunit

        function $hum_op(model::ASHRAEMoistAir,st::ThermodynamicState,unit=$humunit)
            return $hum_op(state_type(st),model,st,unit)
        end

        function $hum_op(mt::HumPT,model::ASHRAEMoistAir,st::ThermodynamicState,unit)
            hum = humidity_spec(FromState(),st)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $hum_op_impl(mt,model,p,t,hum)
            return convert_unit($humunit,unit,res)
        end

        function $hum_op(mt::MultiPT,model::ASHRAEMoistAir,st::ThermodynamicState,unit)
            x = mol_fraction(FromState(),st,nothing,molecular_weight(model))
            hum = last(x)
            _mt = MolHumidityPT() #MolHumidityPT is a complete state type, so it can be called directly instead of using QuickTypes
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $hum_op_impl(_mt,model,p,t,hum)
            return convert_unit($humunit,unit,res)
        end
    end
end


for (op,_unit) in zip([:mol_volume,:mass_volume,
        :mol_density,:mass_density,
        :mol_entropy,:mass_entropy,
        :mol_enthalpy,:mass_enthalpy,
        :compressibility_factor],
        [u"m^3/mol",u"m^3/kg",
        u"mol/(m^3)",u"kg/(m^3)",
        u"J/(mol*K)",u"J/(kg*K)",
        u"J/mol",u"J/kg",
        Unitful.NoUnits])


    op_impl = Symbol(op,:_impl)
    @eval begin
        function $op(model::ASHRAEMoistAir,st::ThermodynamicState,unit=$_unit)
            return $op(state_type(st),model,st,unit)
        end

        function $op(model::ASHRAEWaterVapor,st::ThermodynamicState,unit=$_unit)
            return $op(state_type(st),model,st,unit)
        end

        function $op(model::ASHRAEDryAir,st::ThermodynamicState,unit=$_unit)
            return $op(state_type(st),model,st,unit)
        end

        function $op(mt::SinglePT,model::ASHRAEWaterVapor,st::ThermodynamicState,unit)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $op_impl(mt,model,p,t)
            return convert_unit($_unit,unit,res)
        end

        function $op(mt::SinglePT,model::ASHRAEDryAir,st::ThermodynamicState,unit)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $op_impl(mt,model,p,t)
            return convert_unit($_unit,unit,res)
        end

        function $op(mt::HumPT,model::ASHRAEMoistAir,st::ThermodynamicState,unit)
            hum = humidity_spec(FromState(),st)
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $op_impl(mt,model,p,t,hum)
            return convert_unit($_unit,unit,res)
        end

        function $op(mt::MultiPT,model::ASHRAEMoistAir,st::ThermodynamicState,unit)
            x = mol_fraction(FromState(),st,nothing,molecular_weight(model))
            hum = last(x)
            _mt = MolHumidityPT() #MolHumidityPT is a complete state type, so it can be called directly instead of using QuickTypes
            p = pressure(FromState(),st)
            t = temperature(FromState(),st)
            res  = $op_impl(mt,model,p,t,hum)
            return convert_unit($_unit,unit,res)
        end
    end
end

#total functions
function total_volume(model::ASHRAEMoistAir,st::ThermodynamicState,unit=u"m^3")
    val = mass_volume(model,st)
    m = mass(FromState(),st,u"kg",molecular_weight(model))
    res = val*m
    return convert_unit(u"m^3",unit,res)
end

function total_entropy(model::ASHRAEMoistAir,st::ThermodynamicState,unit=u"J/K")
    val = mass_entropy(model,st)
    m = mass(FromState(),st,u"kg",molecular_weight(model))
    res = val*m
    return convert_unit(u"J/K",unit,res)
end

function total_enthalpy(model:: AbstractMoistAirModel,st::ThermodynamicState,unit=u"J")
    val = mass_enthalpy(model,st)
    m = mass(FromState(),st,u"kg",molecular_weight(model))
    res = val*m
    return convert_unit(u"J",unit,res)
end

function mol_gibbs(model::AbstractMoistAirModel,st::ThermodynamicState,unit=u"J/mol")
    h = mol_enthalpy(model,st)
    t = temperature(FromState(),st)
    s = mol_entropy(model,st)
    res = h-ts
    return convert_unit("J/mol",unit,res)
end

function mass_gibbs(model::AbstractMoistAirModel,st::ThermodynamicState,unit=u"J/kg")
    h = mass_enthalpy(model,st)
    t = temperature(FromState(),st)
    s = mass_entropy(model,st)
    res = h-ts
    return convert_unit("J/kg",unit,res)
end

function total_gibbs(model::AbstractMoistAirModel,st::ThermodynamicState,unit=u"J")
    h = mol_enthalpy(model,st)
    t = temperature(FromState(),st)
    s = mol_entropy(model,st)
    n = moles(FromState,st,u"mol",molecular_weight(model))
    res = n*(h-ts)
    return convert_unit("J",unit,res)
end

function mol_fraction(model::ASHRAEMoistAir,st::ThermodynamicState)
    xv = hum_molfrac(model,st)
    return SVector(1-xv,xv)
     1
end

function mass_fraction(model::ASHRAEMoistAir,st::ThermodynamicState)
    xv = hum_massfrac(model,st)
    return SVector(1-xv,xv)
end

function mol_number(model::ASHRAEMoistAir,st::ThermodynamicState)
    xv = hum_molfrac(model,st)
    n = moles(FromState(),st,u"mol",molecular_weight(model))
    return SVector(n*(1-xv),n*xv)
end

function mass_number(model::ASHRAEMoistAir,st::ThermodynamicState)
    xv = hum_massfrac(model,st)
    n = moles(FromState(),st,u"mol",molecular_weight(model))
    mw = moist_mw(xv)
    m = n*mw
    return SVector(m*(1-xv),m*xv)
end

