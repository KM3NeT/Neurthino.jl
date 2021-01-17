module Neurthino

using LinearAlgebra
using SparseArrays
using StaticArrays
using Polynomials
using DocStringExtensions
using Interpolations
using LRUCache

export transprob, OscillationParameters, PMNSMatrix, Hamiltonian, MatterOscillationMatrices
export masssquareddiff!, cpphase!, mixingangle!

export NeutrinoFlavor, Electron, Muon, Tau

const N_A = 6.022e23 #[mol^-1]
const G_F = 8.961877245622253e-38 #[eV*cm^3]

include("Oscillation.jl")
include("Matter.jl")
include("PREM.jl")

end # module
