module Neurthino

using LinearAlgebra
using SparseArrays

import Polynomials: Poly

include("PREM.jl")

struct OscillationParameters
    mixing_angles::AbstractSparseMatrix{T, S} where {T <: Real, S <: Integer}
    mass_::AbstractSparseMatrix{T, S} where {T <: Real, S <: Integer}
    cp_phase::T where {T <: Real}
end

"""
    transition_probability(U::AbstractArray{T, 2}, H::AbstractVector{S, 1}, L::R)

Calculate the transistion probability between two neutrino flavours

# Arguments
- `U::AbstractArray{T, 2}`: Unitary transistion matrix
- `H::AbstractVector{S}`:   Energy eigenvalues
- `L::R`:                   Baseline

"""
function transition_probability(U::AbstractArray{T, 2}, H::AbstractVector{S}, L::R) where {T <: Number, S <: Real, R <: Real} 
    H_diag = Diagonal(H)
    A = adjoint(U) * exp(-1im * H_diag * L) * U
    P = abs.(A) .^ 2
end


"""
    number_cp_phases(n::Unsigned)

Returns the number of CP violating phases at given number of neutrino types

# Arguments
- `n::Unsigned`: number of neutrino types in the supposed model

# Examples
```julia-repl
julia> Neurthino.number_cp_phases(3)
1
```
"""
function number_cp_phases(n::T) where {T <: Integer}
    if (n < 1) return 0 end
    cp_phases = div( (n-1)*(n-2) , 2 )
end


"""
    number_mixing_angles(n::Unsigned)

Returns the number of mixing angles at given number of neutrino types

# Arguments
- `n::Unsigned`: number of neutrino types in the supposed model

# Examples
```julia-repl
julia> Neurthino.number_mixing_phases(3)
3
```
"""
function number_mixing_angles(n::T) where {T <: Unsigned}
    if (n < 1) return 0 end
    mixing_angles = div( n*(n-1) , 2 )
end



end # module
