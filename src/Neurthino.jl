module Neurthino

using LinearAlgebra
using SparseArrays

import Polynomials: Poly

include("PREM.jl")

struct PMNSMatrix
    matrix::AbstractSparseMatrix{T, S} where {T <: Number, S <: Integer}
end

struct Hamlitionian
    matrix::AbstractSparseMatrix{T, S} where {T <: Number, S <: Integer}
end

struct OscillationParameters
    dim::Integer
    mixing_angles::AbstractSparseMatrix{T, S} where {T <: Real, S <: Integer}
    mass_squared_diff::AbstractSparseMatrix{T, S} where {T <: Real, S <: Integer}
    cp_phases::AbstractSparseMatrix{T, S} where {T <: Real, S <: Integer}
end

function _generate_ordered_index_pairs(n::Integer)
    number_of_angles = number_mixing_angles(n)
    indices = Vector{Pair{Int64,Int64}}(undef, number_of_angles)
    a = 1
    for i in 1:n 
        for j in 1:i-1
            indices[a] = Pair(j,i)
            a += 1
        end
    end
    indices
end

function PMNSMatrix(osc_params::OscillationParameters)
    pmns = sparse(1.0I, osc_params.dim, osc_params.dim) 
    indices = _generate_ordered_index_pairs(osc_params.dim)
    for (i, j) in indices
        rot = sparse(1.0I, osc_params.dim, osc_params.dim) 
        mixing_angle = osc_params.mixing_angles[i, j]
        c, s = cos(mixing_angle), sin(mixing_angle)
        rot[i, i] = c
        rot[j, j] = c
        rot[i, j] = s
        rot[j, i] = -s
        if CartesianIndex(i, j) in osc_params.cp_phases
            cp_phase = osc_params.cp_phases[i, j]
            cp_term = exp(-1im * cp_phase)
            rot[i, j] *= cp_term
            rot[j, i] *= conj(cp_term)
        end
        pmns = rot * pmns 
    end
    PMNSMatrix(pmns)
end

function Hamiltonian(osc_params::OscillationParameters)
    H = spzeros(osc_params.dim, osc_params)
    for (i, j, v) in zip(findnz(osc_params.mass_squared_diff)...)
        if i < j
            H[i, i] += v
        elseif j < i
            H[i, i] -= v
        end
    end
    H /= osc_params.dim
    Hamiltonian(H)
end


"""
    transition_probability(U::AbstractArray{T, 2}, H::AbstractVector{S, 1}, L::R)

Calculate the transistion probability between two neutrino flavours

# Arguments
- `U::AbstractArray{T, 2}`: Unitary transistion matrix
- `H::AbstractVector{S}`:   Energy eigenvalues
- `L::R`:                   Baseline

"""
function transition_probability(U::PMNSMatrix, H::AbstractVector{S}, baseline::R) where {T <: Number, S <: Real, R <: Real} 
    H_diag = Diagonal(H)
    A = adjoint(U) * exp(-1im * H_diag * L) * U
    P = A * adjoint(A) 
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
function number_mixing_angles(n::T) where {T <: Integer}
    if (n < 1) return 0 end
    mixing_angles = div( n*(n-1) , 2 )
end

end # module
