# =================================================================================================
# =                                         Hyland83a.jl                                          =
# =================================================================================================

#=
Testing the virial coefficients. Table 1 from reference [2].
=#
@testset "Hyland83a" begin
@testset "Virial Coeffients" begin

@testset "Baa" begin
@test -55.94 ≈ MoistAir.Baa(173.15)*1e6 atol=0.01
@test -13.15 ≈ MoistAir.Baa(273.15)*1e6 atol=0.01
@test 3.72 ≈ MoistAir.Baa(373.15)*1e6 atol=0.01
@test 12.31 ≈ MoistAir.Baa(473.15)*1e6 atol=0.01
end


@testset "dBaa/dT" begin
@test 0.7240 ≈ MoistAir.dBaa(173.15)*1e6 atol=0.0001
@test 0.2460 ≈ MoistAir.dBaa(273.15)*1e6 atol=0.0001
@test 0.1146 ≈ MoistAir.dBaa(373.15)*1e6 atol=0.0001
@test 0.0640 ≈ MoistAir.dBaa(473.15)*1e6 atol=0.0001
end

@testset "Caaa" begin
@test 2267.0 ≈ MoistAir.Caaa(173.15)*1e12 atol=1.0
@test 1409.0 ≈ MoistAir.Caaa(273.15)*1e12 atol=1.0
@test 1202.0 ≈ MoistAir.Caaa(373.15)*1e12 atol=1.0
@test 1139.0 ≈ MoistAir.Caaa(473.15)*1e12 atol=1.0
end

@testset "dCaaa/dT" begin
@test -18.00 ≈ MoistAir.dCaaa(173.15)*1e12 atol=0.01
@test -3.65 ≈ MoistAir.dCaaa(273.15)*1e12 atol=0.01
@test -1.06 ≈ MoistAir.dCaaa(373.15)*1e12 atol=0.01
@test -0.34 ≈ MoistAir.dCaaa(473.15)*1e12 atol=0.01
end

@testset "Baw" begin
@test -93.3 ≈ MoistAir.Baw(173.15)*1e6 atol=0.1
@test -36.4 ≈ MoistAir.Baw(273.15)*1e6 atol=0.1
@test -14.5 ≈ MoistAir.Baw(373.15)*1e6 atol=0.1
@test -3.1 ≈ MoistAir.Baw(473.15)*1e6 atol=0.1
end

@testset "dBaw/dT" begin
@test 1.01 ≈ MoistAir.dBaw(173.15)*1e6 atol=0.01
@test 0.318 ≈ MoistAir.dBaw(273.15)*1e6 atol=0.001
@test 0.151 ≈ MoistAir.dBaw(373.15)*1e6 atol=0.001
@test 0.0869 ≈ MoistAir.dBaw(473.15)*1e6 atol=0.0001
end

@testset "Caaw" begin
@test 1023 ≈ MoistAir.Caaw(173.15)*1e12 atol=1.0
@test 861 ≈ MoistAir.Caaw(273.15)*1e12 atol=1.0
@test 696 ≈ MoistAir.Caaw(373.15)*1e12 atol=1.0
@test 627 ≈ MoistAir.Caaw(473.15)*1e12 atol=1.0
end

@testset "dCaaw/dT" begin
@test  5.56 ≈ MoistAir.dCaaw(173.15)*1e12 atol=0.01
@test -2.44 ≈ MoistAir.dCaaw(273.15)*1e12 atol=0.01
@test -1.02 ≈ MoistAir.dCaaw(373.15)*1e12 atol=0.01
@test -0.457 ≈ MoistAir.dCaaw(473.15)*1e12 atol=0.001
end

@testset "Caww" begin
@test -2.0e7 ≈ MoistAir.Caww(173.15)*1e12 atol=0.1e7
@test -2.2e5 ≈ MoistAir.Caww(273.15)*1e12 atol=0.1e5
@test -3.0e4 ≈ MoistAir.Caww(373.15)*1e12 atol=0.1e4
@test -8.4e3 ≈ MoistAir.Caww(473.15)*1e12 atol=0.1e3
end

@testset "dCaww/dT" begin
@test 1.61e6 ≈ MoistAir.dCaww(173.15)*1e12 atol=0.01e6
@test 6.05e3 ≈ MoistAir.dCaww(273.15)*1e12 atol=0.01e3
@test    456 ≈ MoistAir.dCaww(373.15)*1e12 atol=1.0
@test   86.9 ≈ MoistAir.dCaww(473.15)*1e12 atol=0.1
# Note that on table 1, the last line is 8.69 instead of 86.9. Probably a typo in the paper.
end
end

@testset "Specific volume of dry air" begin
# Table 3
P = 0.1e6
@test 0.49510 ≈ MoistAir.volumeair(173.15, P) atol=0.00003
@test 0.78364 ≈ MoistAir.volumeair(273.15, P) atol=0.0001
@test 1.0713 ≈ MoistAir.volumeair(373.15, P) atol=0.0001
@test 1.3586 ≈ MoistAir.volumeair(473.15, P) atol=0.0001

P = 0.5e6
@test 0.09747 ≈ MoistAir.volumeair(173.15, P) atol=0.00003
@test 0.15637 ≈ MoistAir.volumeair(273.15, P) atol=0.00003
@test 0.21436 ≈ MoistAir.volumeair(373.15, P) atol=0.00003
@test 0.27207 ≈ MoistAir.volumeair(473.15, P) atol=0.00003

P = 1e6
@test 0.04775 ≈ MoistAir.volumeair(173.15, P) atol=0.00003
@test 0.07797 ≈ MoistAir.volumeair(273.15, P) atol=0.00003
@test 0.10725 ≈ MoistAir.volumeair(373.15, P) atol=0.00003
@test 0.13625 ≈ MoistAir.volumeair(473.15, P) atol=0.00003

P = 5e6
@test 0.00795 ≈ MoistAir.volumeair(173.15, P) atol=0.00003
@test 0.01533 ≈ MoistAir.volumeair(273.15, P) atol=0.00003
@test 0.02162 ≈ MoistAir.volumeair(373.15, P) atol=0.00003
@test 0.02763 ≈ MoistAir.volumeair(473.15, P) atol=0.00003
end

@testset "Specific enthalpy of dry air" begin
# Table 3
P = 0.1e6
@test -100.63 ≈ MoistAir.enthalpyair(173.15, P)/1000 atol=0.01
@test    0.004≈ MoistAir.enthalpyair(273.15, P)/1000 atol=0.01
@test  100.79 ≈ MoistAir.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.55 ≈ MoistAir.enthalpyair(473.15, P)/1000 atol=0.01

P = 0.5e6
@test -103.17 ≈ MoistAir.enthalpyair(173.15, P)/1000 atol=0.01
@test   -1.10 ≈ MoistAir.enthalpyair(273.15, P)/1000 atol=0.01
@test  100.26 ≈ MoistAir.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.30 ≈ MoistAir.enthalpyair(473.15, P)/1000 atol=0.01

P = 1e6
@test -106.42 ≈ MoistAir.enthalpyair(173.15, P)/1000 atol=0.01
@test   -2.48 ≈ MoistAir.enthalpyair(273.15, P)/1000 atol=0.01
@test   99.60 ≈ MoistAir.enthalpyair(373.15, P)/1000 atol=0.01
@test  202.00 ≈ MoistAir.enthalpyair(473.15, P)/1000 atol=0.01

P = 5e6
@test -135.55 ≈ MoistAir.enthalpyair(173.15, P)/1000 atol=0.01
@test  -13.15 ≈ MoistAir.enthalpyair(273.15, P)/1000 atol=0.01
@test   94.63 ≈ MoistAir.enthalpyair(373.15, P)/1000 atol=0.01
@test  199.8 ≈ MoistAir.enthalpyair(473.15, P)/1000 atol=0.1
end


@testset "Specific entropy of dry air" begin
P = 0.1e6
@test -0.4550 ≈ MoistAir.entropyair(173.15, P)/1000 atol=0.0001
@test  0.0038 ≈ MoistAir.entropyair(273.15, P)/1000 atol=0.0001
@test  0.3182 ≈ MoistAir.entropyair(373.15, P)/1000 atol=0.0001
@test  0.5597 ≈ MoistAir.entropyair(473.15, P)/1000 atol=0.0001
end

# Table 4, appendix of reference [2]
@testset "Volume of saturated moist air" begin
volmoist(Tk, P) = MoistAir.volumemoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)
P = 0.1e6
@test 0.49510 ≈ volmoist(173.15, P) atol=0.00005
@test 0.78846 ≈ volmoist(273.15, P) atol=0.00005
@test 3.4974 ≈ volmoist(363.15, P) atol=0.0002 #fail
P = 0.5e6
@test 0.097466 ≈ volmoist(173.15, P) atol=0.000006
@test 0.15657 ≈ volmoist(273.15, P) atol=0.00005
@test 0.24273 ≈ volmoist(363.15, P) atol=0.00003
@test 0.26904 ≈ volmoist(373.15, P) atol=0.00003
P=1e6
@test 0.047752 ≈ volmoist(173.15, P) atol=0.000006
@test 0.078022 ≈ volmoist(273.15, P) atol=0.000005
@test 0.11226 ≈ volmoist(363.15, P) atol=0.00001
@test 0.11942 ≈ volmoist(373.15, P) atol=0.00001
@test 0.23281 ≈ volmoist(423.15, P) atol=0.00002 #fail
@test 0.32711 ≈ volmoist(433.15, P) atol=0.00002
@test 0.61489 ≈ volmoist(443.15, P) atol=0.00005 #fail

P=5e6
@test 0.007951 ≈ volmoist(173.15, P) atol=0.000006
@test 0.015331 ≈ volmoist(273.15, P) atol=0.000005
@test 0.021310 ≈ volmoist(363.15, P) atol=0.00001
@test 0.022073 ≈ volmoist(373.15, P) atol=0.00001
@test 0.02727  ≈ volmoist(423.15, P) atol=0.00001
@test 0.02885  ≈ volmoist(433.15, P) atol=0.00001
@test 0.03077  ≈ volmoist(443.15, P) atol=0.00001
@test 0.03315  ≈ volmoist(453.15, P) atol=0.00001
@test 0.0402  ≈ volmoist(473.15, P) atol=0.0001
end

@testset "enthalpy of saturated moist air" begin
hmoist(Tk, P) = MoistAir.enthalpymoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)/1000
P = 0.1e6
@test -100.627 ≈ hmoist(173.15, P) atol=0.005
@test    9.602 ≈ hmoist(273.15, P) atol=0.005
@test   4035.0 ≈ hmoist(363.15, P) atol=0.3 # In the paper the number is 4305! #fail

P = 0.5e6
@test -103.168 ≈ hmoist(173.15, P) atol=0.005
@test    0.831 ≈ hmoist(273.15, P) atol=0.005
@test   365.55 ≈ hmoist(363.15, P) atol=0.03
@test   533.27 ≈ hmoist(373.15, P) atol=0.03
@test   35957 ≈ hmoist(423.15, P) atol=2.0 #fail

P = 1e6
@test -106.415 ≈ hmoist(173.15, P) atol=0.005
@test   -1.498 ≈ hmoist(273.15, P) atol=0.005
@test   217.77 ≈ hmoist(363.15, P) atol=0.03
@test   293.17 ≈ hmoist(373.15, P) atol=0.03
@test   1788.8 ≈ hmoist(423.15, P) atol=0.5 #fail
@test   3113.3 ≈ hmoist(433.15, P) atol=0.5
@test   7206.0 ≈ hmoist(443.15, P) atol=2.0 #fail

P = 5e6
@test -135.547 ≈ hmoist(173.15, P) atol=0.005
@test  -12.929 ≈ hmoist(273.15, P) atol=0.005
@test   109.89 ≈ hmoist(363.15, P) atol=0.03
@test   132.42 ≈ hmoist(373.15, P) atol=0.03
@test   347.3 ≈ hmoist(423.15, P) atol=0.5
@test   428.8 ≈ hmoist(433.15, P) atol=0.5
@test   534.1 ≈ hmoist(443.15, P) atol=0.5
@test   672.0 ≈ hmoist(453.15, P) atol=1.0
@test   1113.0 ≈ hmoist(473.15, P) atol=1.0
end

@testset "Entropy of saturated moist air" begin
# Table 4 in the appendix of reference [2].
smoist(Tk, P) = MoistAir.entropymoist(Tk, P, MoistAir.efactor(Tk,P)*MoistAir.Pws(Tk)/P)/1000
P = 0.1e6
@test  -0.45499 ≈ smoist(173.15, P) atol=0.00003
@test   0.04070 ≈ smoist(273.15, P) atol=0.00003
@test  11.7287 ≈ smoist(363.15, P) atol=0.0005 #fail

P = 0.5e6
@test  -0.92718 ≈ smoist(173.15, P) atol=0.00005
@test  -0.45417 ≈ smoist(273.15, P) atol=0.00005
@test   0.64517 ≈ smoist(363.15, P) atol=0.0005
@test   1.1102 ≈ smoist(373.15, P) atol=0.0005
@test  90.01 ≈ smoist(423.15, P) atol=0.05 #fail

P = 1e6
@test  -1.13929 ≈ smoist(173.15, P) atol=0.00005 # In the table the signal is +
@test  -0.66104 ≈ smoist(273.15, P) atol=0.00005
@test   0.00805 ≈ smoist(363.15, P) atol=0.00005
@test   0.21676 ≈ smoist(373.15, P) atol=0.0005
@test   4.0478  ≈ smoist(423.15, P) atol=0.0005 #fail
@test   7.3103  ≈ smoist(433.15, P) atol=0.0005
@test   17.221  ≈ smoist(443.15, P) atol=0.005 #fail

P = 5e6
@test  -1.72380 ≈ smoist(173.15, P) atol=0.0001 # In the table the signal is +
@test  -1.15924 ≈ smoist(273.15, P) atol=0.00005
@test  -0.77519 ≈ smoist(363.15, P) atol=0.00005
@test  -0.71331 ≈ smoist(373.15, P) atol=0.00005
@test  -0.1629  ≈ smoist(423.15, P) atol=0.0002
@test   0.0364  ≈ smoist(433.15, P) atol=0.0002
@test   0.2901  ≈ smoist(443.15, P) atol=0.0002
@test   0.6185  ≈ smoist(453.15, P) atol=0.0002
@test   1.642  ≈ smoist(473.15, P) atol=0.002
end
end