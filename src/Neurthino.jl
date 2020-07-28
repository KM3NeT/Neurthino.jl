module Neurthino

using LinearAlgebra
using SparseArrays
using PhysicalConstants
using Unitful

using Polynomials

# TEMPORARY UNTIL PhysicalConstants.jl GETS UPDATED
G_F = 1.1663787e-5u"GeV^-2"
G_F = G_F * (PhysicalConstants.CODATA2018.SpeedOfLightInVacuum * PhysicalConstants.CODATA2018.ReducedPlanckConstant)^3
G_F = uconvert(u"eV*cm^3", G_F)


include("PREM.jl")
include("Oscillation.jl")

end # module
