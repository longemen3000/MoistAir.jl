Tk = 183.3:10.0:463.3
Pa = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 49.0]*1e5
xvs(Tk, P) = MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P

volmoist(Tk, P) = MoistAir.volumemoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)

@testset "volume - saturated moist air" begin
    st1 = state(t=λ,p=λ,rel_hum=1.0)
    for p1 in Pa, t1 in Tk
        global x = xvs(t1, p1)
        if x <= 1.0
            v1 = volmoist(t1, p1)
            v2 = dry_volume(MoistAirModel(), st1(t1,p1))
            @test v1 ≈ v2 rtol=1e-6
        end
    end
end

@testset "volume - dry air" begin
    st1 = state(t=λ,p=λ,rel_hum=0.0)
    st2 = state(t=λ,p=λ)
    for p1 in Pa, t1 in Tk
        v1 = MoistAir.volumeair(t1, p1)
        v2 = MoistAir.dry_volume(MoistAirModel(),st1(t1,p1))
        v3 = MoistAir.mass_volume(ASHRAEDryAir(), st2(t1, p1))
        @test v1 ≈ v2 rtol=1e-6
        @test v2 ≈ v3 rtol=1e-6
    end
end

@testset "enthalpy - saturated moist air" begin
    hmoist(Tk, P) = MoistAir.enthalpymoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)/1000
    st1 = state(t=λ,p=λ,rel_hum=1.0)
    for p1 in Pa, t1 in Tk
        global x = xvs(t1, p1)
        if x <= 1.0
            v1 = hmoist(t1, p1)
            v2 = dry_enthalpy(MoistAirModel(), st1(t1,p1))*0.001
            @test v1 ≈ v2 rtol=1e-4
        end
    end
end

@testset "enthalpy - dry air" begin
    st1 = state(t=λ,p=λ,rel_hum=0.0)
    st2 = state(t=λ,p=λ)
    for p1 in Pa, t1 in Tk
        v1 = MoistAir.enthalpyair(t1, p1)
        v2 = MoistAir.dry_enthalpy(MoistAirModel(),st1(t1,p1))
        v3 = MoistAir.mass_enthalpy(ASHRAEDryAir(), st2(t1, p1))    #println(p1, " - ", t1, " - ", v1, " - ", v2, " - ", v3)
        @test v1 ≈ v2 rtol=1e-4
        @test v2 ≈ v3 rtol=1e-4
    end
end

@testset "entropy - saturated moist air" begin
smoist(Tk, P) = MoistAir.entropymoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)/1000
st1 = state(t=λ,p=λ,rel_hum=1.0)
for p1 in Pa, t1 in Tk
    global x = xvs(t1, p1)
    if x <= 1.0
        v1 = smoist(t1, p1)
        v2 = dry_entropy(MoistAirModel(), st1(t1,p1))*0.001
        @test v1 ≈ v2 rtol=1e-6
    end
end
end

@testset "entropy - dry air" begin
st1 = state(t=λ,p=λ,rel_hum=0.0)
st2 = state(t=λ,p=λ)
for p1 in Pa, t1 in Tk
    v1 = MoistAir.entropyair(t1, p1)
    v2 = MoistAir.dry_entropy(MoistAirModel(),st1(t1,p1))
    v3 = MoistAir.mass_entropy(ASHRAEDryAir(), st2(t1, p1)) 
    #println(p1, " - ", t1, " - ", v1, " - ", v2, " - ", v3)
    @test v1 ≈ v2 rtol=1e-4
    @test v2 ≈ v3 rtol=1e-4
end
end

@testset "Calculations of M* in table 13 of reference [5]" begin
    function mfun(Tbu, P=101325.0)
        st1 = state(t=λ,p=λ,rel_hum=1.0)
        hs = MoistAir.dry_enthalpy(MoistAirModel(),st1(Tbu,P))
        ws = MoistAir.hum_ratio(MoistAirModel(), st1(Tbu,P))
        hw = MoistAir.enthalpywi(Tbu)
        return (hs - ws*hw)/1000
    end


    @test 12.94 ≈ mfun(2.0+273.15) atol=0.015
    @test 20.50 ≈ mfun(6.0+273.15) atol=0.015
    @test 29.03 ≈ mfun(10.0+273.15) atol=0.015
    @test 38.78 ≈ mfun(14.0+273.15) atol=0.015
    @test 50.03 ≈ mfun(18.0+273.15) atol=0.015
    @test 63.11 ≈ mfun(22.0+273.15) atol=0.015
    @test 78.46 ≈ mfun(26.0+273.15) atol=0.015
    @test 96.57 ≈ mfun(30.0+273.15) atol=0.015
    @test 118.07 ≈ mfun(34.0+273.15) atol=0.015
    @test 143.74 ≈ mfun(38.0+273.15) atol=0.015
    @test 174.58 ≈ mfun(42.0+273.15) atol=0.015
    @test 192.31 ≈ mfun(44.0+273.15) atol=0.015
    @test 233.36 ≈ mfun(48.0+273.15) atol=0.015

end

@testset "Testing Table 20 from reference [5]" begin

    Tk = vcat(5*ones(4), 25*ones(4), 50*ones(4)) .+ 273.15
    Tbu = [5.0, 2.0, -1.0, -3.0, 25.0, 20.0, 15.0, 10.0, 25.0, 22.0, 20.0, 19.0] .+ 273.15
    Td = [5.0, -2.16, -11.92, -37.23, 25.0, 17.60, 7.73, -10.42, 13.47, 4.16, -5.48, -14.12]
    w = [5.42, 3.16, 1.35, 0.11, 20.17, 12.66, 6.56, 1.55, 9.67, 5.11, 2.39, 1.11]*1e-3
    rel = [100.0, 58.6, 25.1, 2.0, 100, 63.5, 33.2, 7.9, 12.5, 6.7, 3.1, 1.5]
    vol = [0.794, 0.792, 0.789, 0.788, 0.872, 0.862, 0.853, 0.846, 0.930, 0.923, 0.919, 0.917]
    h = [18.64, 12.97, 8.42, 5.30, 76.50, 57.38, 41.85, 29.09, 75.40, 63.58, 56.51, 53.20]
    s = [0.0697, 0.0490, 0.0319, 0.0195, 0.2698, 0.2048, 0.1506, 0.1039, 0.2610, 0.2192, 0.1934, 0.1808]

    stw = state(t=λ,p=λ,hum_wetbulb=λ)
    P = 101325.0
    Td1 = MoistAir.hum_dewpoint.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu)).- 273.15
    w1 =  MoistAir.hum_ratio.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu))
    rel1 =  MoistAir.rel_hum.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu)) * 100
    vol1 =  MoistAir.dry_volume.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu)) 
    h1 =  MoistAir.dry_enthalpy.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu)) /1000
    s1 =  MoistAir.dry_entropy.(Ref(MoistAirModel()), stw.(Tk,Ref(P),Tbu)) /1000

    @test maximum(abs, Td1-Td) ≈ 0.0 atol=0.01
    @test maximum(abs, w1-w) ≈ 0.0 atol=0.00001
    @test maximum(abs, rel1-rel) ≈ 0.0 atol=0.1
    @test maximum(abs, vol1-vol) ≈ 0.0 atol=0.001
    @test maximum(abs, h1-h) ≈ 0.0 atol=0.01
    @test maximum(abs, s1-s) ≈ 0.0 atol=0.0001

    stbu1 = state(t=λ,p=λ,hum_dewpoint=λ)
    Tbu1 = MoistAir.hum_wetbulb.(Ref(MoistAirModel()), stbu1.(Tk,P,Td1 .+ 273.15))
    @test maximum(abs, Tbu1-Tbu) ≈ 0.0 atol=1e-4

    stbu2 = state(t=λ,p=λ,rel_hum=λ) 
    Tbu2 = MoistAir.hum_wetbulb.(Ref(MoistAirModel()), stbu2.(Tk,P,rel1/100))
    @test maximum(abs, Tbu2-Tbu) ≈ 0.0 atol=1e-4

    stbu3 = state(t=λ,p=λ,hum_ratio=λ) 
    Tbu3 = MoistAir.hum_wetbulb.(Ref(MoistAirModel()),stbu3.(Tk,P,w1/100))
    @test maximum(abs, Tbu2-Tbu) ≈ 0.0 atol=1e-4


    std2 = state(t=λ,p=λ,hum_ratio=λ) 
    Td2 = MoistAir.hum_dewpoint.(Ref(MoistAirModel()), std2.(Tk,P,w1)) .- 273.15
    @test maximum(abs, Td2-Td1) ≈ 0.0 atol=1e-4


    std3 = state(t=λ,p=λ,rel_hum=λ) 
    Td3 = MoistAir.hum_dewpoint.(Ref(MoistAirModel()), std3.(Tk,P,rel1/100)) .- 273.15
    @test maximum(abs, Td3-Td1) ≈ 0.0 atol=1e-4

end
