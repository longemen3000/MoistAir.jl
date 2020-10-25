# =================================================================================================
# =                                         Hyland83.jl                                           =
# =================================================================================================
@testset "Hyland83" begin
# Appendix of ref. [1]. First table
@testset "Specific volume of saturated ice" begin
@test 1.0768e-3 ≈ MoistAir.volumeice(173.15) atol=0.0001e-3
@test 1.0829e-3 ≈ MoistAir.volumeice(223.15) atol=0.0001e-3
@test 1.0909e-3 ≈ MoistAir.volumeice(273.16) atol=0.0001e-3
end
@testset "Specific volume of saturated liquid water" begin
@test 1.00021e-3 ≈ MoistAir.volumewater(273.16) atol=0.00001e-3
@test 1.01215e-3 ≈ MoistAir.volumewater(323.15) atol=0.00001e-3
@test 1.04346e-3 ≈ MoistAir.volumewater(373.15) atol=0.00001e-3 # In the paper there is a typo in the power of 10: 1.04346e10
@test 1.09050e-3 ≈ MoistAir.volumewater(423.15) atol=0.00001e-3
@test 1.15653e-3 ≈ MoistAir.volumewater(473.15) atol=0.00001e-3
end

@testset "Saturation enthalpy of saturated ice" begin
@test -507.215e3 ≈ MoistAir.enthalpyice(173.15) atol=0.001e3
@test -429.413e3 ≈ MoistAir.enthalpyice(223.15) atol=0.001e3
@test -333.429e3 ≈ MoistAir.enthalpyice(273.15) atol=0.001e3  # The table presents -333.409. Typo???
end

@testset "Saturation enthalpy of saturated liquid water" begin
@test 0.0 ≈ MoistAir.enthalpywater(273.16) atol=1.0
@test 209.330e3 ≈ MoistAir.enthalpywater(323.15) atol=5
@test 419.158e3 ≈ MoistAir.enthalpywater(373.15) atol=5
@test 632.210e3 ≈ MoistAir.enthalpywater(423.15) atol=9
@test 852.329e3 ≈ MoistAir.enthalpywater(473.15) atol=15
end

# Table 2 of reference [2]
@testset "Enhancement factor" begin
P = 0.1e6
@test 1.0105 ≈ MoistAir.efactor(173.15, P) atol=0.0001
@test 1.0039 ≈ MoistAir.efactor(273.15, P) atol=0.0001
@test 1.0039 ≈ MoistAir.efactor(363.15, P) atol=0.0001
P = 0.5e6
@test 1.054 ≈ MoistAir.efactor(173.15, P) atol=0.001
@test 1.0177 ≈ MoistAir.efactor(273.15, P) atol=0.0006  # I think the table has a typo!!?
@test 1.0180 ≈ MoistAir.efactor(363.15, P) atol=0.0001
@test 1.0188 ≈ MoistAir.efactor(373.15, P) atol=0.0001
@test 1.0022 ≈ MoistAir.efactor(423.15, P) atol=0.0001
P = 1.0e6
@test 1.113 ≈ MoistAir.efactor(173.15, P) atol=0.001
@test 1.0353 ≈ MoistAir.efactor(273.15, P) atol=0.0015  # I think the table has a typo!!?
@test 1.0284 ≈ MoistAir.efactor(363.15, P) atol=0.0001
@test 1.0295 ≈ MoistAir.efactor(373.15, P) atol=0.0001
@test 1.0288 ≈ MoistAir.efactor(423.15, P) atol=0.0001
@test 1.0235 ≈ MoistAir.efactor(433.15, P) atol=0.0001
@test 1.0142 ≈ MoistAir.efactor(443.15, P) atol=0.0001


P = 5.0e6
@test 1.820 ≈ MoistAir.efactor(173.15, P) atol=0.001
@test 1.191 ≈ MoistAir.efactor(273.15, P) atol=0.008  # I think the table has a typo!!?
@test 1.1102 ≈ MoistAir.efactor(363.15, P) atol=0.0001
@test 1.1082 ≈ MoistAir.efactor(373.15, P) atol=0.0001
@test 1.111 ≈ MoistAir.efactor(423.15, P) atol=0.001
@test 1.114 ≈ MoistAir.efactor(433.15, P) atol=0.001
@test 1.116 ≈ MoistAir.efactor(443.15, P) atol=0.001
@test 1.117 ≈ MoistAir.efactor(453.15, P) atol=0.001
@test 1.116 ≈ MoistAir.efactor(473.15, P) atol=0.001
end

@testset "Saturation pressure of water vapor over ice" begin
# Second table in the appendix of reference [1]
@test 1.40510e-3 ≈ MoistAir.Pws_s(173.15) atol=0.00001e-3
@test 6.11153e2 ≈ MoistAir.Pws_s(273.15) atol=0.00001e2
end

@testset "Saturation pressure of water vapor over liquid water" begin
# Second table in the appendix of reference [1]
@test 6.11213e2 ≈ MoistAir.Pws_l(273.15) atol=0.00001e2
@test 1.01419e5 ≈ MoistAir.Pws_l(373.15) atol=0.00001e5
@test 1.55507e6 ≈ MoistAir.Pws_l(473.15) atol=0.00001e6
end


@testset " saturation temperature function" begin
@test MoistAir.Tws(MoistAir.Pws_s(173.15)) ≈ 173.15 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_s(223.15)) ≈ 223.15 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_s(273.15)) ≈ 273.15 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_s(273.16)) ≈ 273.16 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_l(273.16)) ≈ 273.16 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_l(323.15)) ≈ 323.15 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_l(373.15)) ≈ 373.15 atol=1e-5
@test MoistAir.Tws(MoistAir.Pws_l(473.15)) ≈ 473.15 atol=1e-5
end

@testset "First virial coefficient of saturated vapor B'" begin
# Second table in the appendix of reference [1]
@test -3.2939e-5 ≈ MoistAir.Blin(173.15) atol=0.0001e-5
@test -8.3497e-7 ≈ MoistAir.Blin(273.15) atol=0.0001e-7
@test -1.4658e-7 ≈ MoistAir.Blin(373.15) atol=0.0001e-5
@test -5.0508e-8 ≈ MoistAir.Blin(473.15) atol=0.0001e-5
end

@testset "Second virial coefficient of saturated vapor C'" begin
# Second table in the appendix of reference [1]
@test -4.6563e-9 ≈ MoistAir.Clin(173.15) atol=0.0001e-9
@test -2.0928e-12 ≈ MoistAir.Clin(273.15) atol=0.0001e-12
@test -5.7548e-14 ≈ MoistAir.Clin(373.15) atol=0.0001e-14
@test -6.3933e-15 ≈ MoistAir.Clin(473.15) atol=0.0001e-15
end

@testset "Specific enthalpy of saturated water vapor" begin
# First table of the appendix of ref. [1]
@test 2315.87 ≈ MoistAir.enthalpyvapor(173.15)/1000 atol=0.01
@test 2408.41 ≈ MoistAir.enthalpyvapor(223.15)/1000 atol=0.01
@test 2500.81 ≈ MoistAir.enthalpyvapor(273.16)/1000 atol=0.01
@test 2591.29 ≈ MoistAir.enthalpyvapor(323.15)/1000 atol=0.01
@test 2675.46 ≈ MoistAir.enthalpyvapor(373.15)/1000 atol=0.01
@test 2746.15 ≈ MoistAir.enthalpyvapor(423.15)/1000 atol=0.01
@test 2793.11 ≈ MoistAir.enthalpyvapor(473.15)/1000 atol=0.01
end

@testset "Specific volume of saturated water vapor" begin
# First table of the appendix of ref. [1]
@test 5.6873e7 ≈ MoistAir.volumevapor(173.15) atol=0.0001e7
@test 2.6146e4 ≈ MoistAir.volumevapor(223.15) atol=0.0001e4
@test 2.0601e2 ≈ MoistAir.volumevapor(273.16) atol=0.0001e2
@test 1.2030e1 ≈ MoistAir.volumevapor(323.15) atol=0.0001e1
@test 1.6718e0 ≈ MoistAir.volumevapor(373.15) atol=0.0001e0
@test 3.9253e-1 ≈ MoistAir.volumevapor(423.15) atol=0.0001e-1
@test 1.2722e-1 ≈ MoistAir.volumevapor(473.15) atol=0.0001e0
end

@testset "Specific entropy of saturated water vapor" begin
# First table of the appendix of ref. [1]
@test 14.30387 ≈ MoistAir.entropyvapor(173.15)/1000 atol=0.00005
@test 11.10965 ≈ MoistAir.entropyvapor(223.15)/1000 atol=0.00005
@test  9.15510 ≈ MoistAir.entropyvapor(273.16)/1000 atol=0.00005
@test  8.07477 ≈ MoistAir.entropyvapor(323.15)/1000 atol=0.00005
@test  7.35365 ≈ MoistAir.entropyvapor(373.15)/1000 atol=0.00005
@test  6.83731 ≈ MoistAir.entropyvapor(423.15)/1000 atol=0.00005
@test  6.43218 ≈ MoistAir.entropyvapor(473.15)/1000 atol=0.00005
end
end