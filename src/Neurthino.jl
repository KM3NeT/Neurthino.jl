module Neurthino

using LinearAlgebra
using SparseArrays
using StaticArrays
using Polynomials
using DocStringExtensions
using Interpolations

export transprob, OscillationParameters, PMNSMatrix, Hamiltonian, MatterOscillationMatrices
export masssquareddiff!, cpphase!, mixingangle!

const N_A = 6.022e23 #[mol^-1]
const G_F = 8.961877245622253e-38 #[eV*cm^3]

include("PREM.jl")
include("Oscillation.jl")
include("Path.jl")
include("Optimization.jl")

end # module
