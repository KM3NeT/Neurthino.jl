using Test

@testset "PREM" begin
include("PREM_tests.jl")
end

@testset "Oscillation" begin
include("oscillation_tests.jl")
end

@testset "Path" begin
include("path_tests.jl")
end

@testset "Matter" begin
include("matter_tests.jl")
end

@testset "Examples" begin
include("examples_tests.jl")
end