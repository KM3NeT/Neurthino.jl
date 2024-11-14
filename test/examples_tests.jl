# check if examples run without error or warning

using Neurthino
using Unitful
using Plots
using LaTeXStrings

include("../examples/matter_oscillogram.jl")
@test_nowarn example_matter_oscillogram()