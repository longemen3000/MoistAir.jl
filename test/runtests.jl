using MoistAir, Unitful,ThermoState
using Test
const λ = VariableSpec()

include("test_hyland83.jl")
include("test_hyland83a.jl")
include("test_user.jl")

