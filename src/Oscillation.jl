@enum NeutrinoFlavor begin
  Electron = 1
  Muon = 2
  Tau = 3
end


struct OscillationParameters{T}
    mixing_angles::Array{T,2}
    mass_squared_diff::Array{T,2}
    cp_phases::Array{T,2}
    dim::Int64

    OscillationParameters(dim::Int64) = begin
        new{Float64}(
                zeros(dim, dim),
                zeros(dim, dim),
                zeros(dim, dim),
                dim)
    end
end

function _generate_ordered_index_pairs(n::Integer)
    indices = Vector{Pair{Int64,Int64}}(undef, mixingangles(n))
    a = 1
    for i in 1:n 
        for j in 1:i-1
            indices[a] = Pair(j,i)
            a += 1
        end
    end
    indices
end

"""
$(SIGNATURES)

Set a mixing angle of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Pair{<:Integer, <:Integer}`: The indices of the mixing angle
- `value` The value which should be applied to the oscillation parameters

"""
function mixingangle!(osc::OscillationParameters, indices::Pair{T, T}, value) where {T <: Integer}
    fromidx = first(indices)
    toidx = last(indices)
    if fromidx < toidx
        osc.mixing_angles[fromidx, toidx] = value
    else
        osc.mixing_angles[toidx, fromidx] = value
    end
end

"""
$(SIGNATURES)

Set a mass squared difference of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Pair{<:Integer, <:Integer}`: The indices of the mass squared difference
- `value` The value which should be applied to the oscillation parameters

"""
function masssquareddiff!(osc::OscillationParameters, indices::Pair{T, T}, value) where {T <: Integer}
    fromidx = first(indices)
    toidx = last(indices)
    if fromidx < toidx
        osc.mass_squared_diff[fromidx, toidx] = value
    else
        osc.mass_squared_diff[toidx, fromidx] = -value
    end
end

"""
$(SIGNATURES)

Set a CP phase of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Pair{<:Integer, <:Integer}`: The indices of the mass difference
- `value` The value which should be applied to the oscillation parameters

"""
function cpphase!(osc::OscillationParameters, indices::Pair{T, T}, value) where {T <: Integer}
    fromidx = first(indices)
    toidx = last(indices)
    if fromidx < toidx
        osc.cp_phases[fromidx, toidx] = value
    else
        osc.cp_phases[toidx, fromidx] = value
    end
end



function PMNSMatrix(osc_params::OscillationParameters)
"""
$(SIGNATURES)

Create rotation matrix (PMNS) based on the given oscillation parameters

# Arguments
- `osc_params::OscillationParameters`: Oscillation parameters

"""
    dim = size(osc_params.mixing_angles)[1]
    pmns = Matrix{ComplexF64}(1.0I, dim, dim) 
    indices = _generate_ordered_index_pairs(dim)
    for (i, j) in indices
        rot = sparse((1.0+0im)I, dim, dim) 
        mixing_angle = osc_params.mixing_angles[i, j]
        c, s = cos(mixing_angle), sin(mixing_angle)
        rot[i, i] = c
        rot[j, j] = c
        rot[i, j] = s
        rot[j, i] = -s
        if CartesianIndex(i, j) in findall(!iszero, osc_params.cp_phases)
            cp_phase = osc_params.cp_phases[i, j]
            cp_term = exp(-1im * cp_phase)
            rot[i, j] *= cp_term
            rot[j, i] *= conj(cp_term)
        end
        pmns = rot * pmns 
    end
    pmns
end

function Hamiltonian(osc_params::OscillationParameters)
"""
$(SIGNATURES)

Create modified hamiltonian matrix consisting of the squared mass differences
based on the given oscillation parameters

# Arguments
- `osc_params::OscillationParameters`: Oscillation parameters

"""
    Hamiltonian(osc_params, zeros(Float64, osc_params.dim))
end 

function Hamiltonian(osc_params::OscillationParameters, lambda)
"""
$(SIGNATURES)

Create modified hamiltonian matrix consisting of the squared mass differences
based on the given oscillation parameters

# Arguments
- `osc_params::OscillationParameters`:  Oscillation parameters
- `lambda`:                             Decay parameters for each mass eigenstate

"""
    H = zeros(ComplexF64, osc_params.dim)
    for i in 1:osc_params.dim
        for j in 1:osc_params.dim
            if i < j
                H[i] += osc_params.mass_squared_diff[i,j]
            elseif j < i
                H[i] -= osc_params.mass_squared_diff[j,i]
            end
        end
        H[i] += 1im * lambda[i]
    end
    H /= osc_params.dim
    H
end


function _transprobampl(U, H, energy, baseline)  
    H_diag = 2.534 * Diagonal{ComplexF64}(H) * baseline / energy 
    U * exp(-1im * H_diag) * adjoint(U)
end

"""
$(SIGNATURES)

Calculate the transistion probabilities between the neutrino flavours

# Arguments
- `U`:          Unitary transition matrix
- `H`:          Energy eigenvalues
- `energy`:     Baseline [km]
- `baseline`:   Energy [GeV]

"""
function transprob(U, H, energy, baseline)  
    A = _transprobampl(U, H, energy, baseline)
    P = abs.(A) .^ 2
end

"""
$(SIGNATURES)

Calculate the transistion probabilities between the neutrino flavours

# Arguments
- `osc_params::OscillationParameters`:  Oscillation parameters
- `energy`:                             Baseline [km]
- `baseline`:                           Energy [GeV]

"""
function transprob(osc_params::OscillationParameters, energy, baseline)  
    H = Hamiltonian(osc_params)
    U = PMNSMatrix(osc_params)
    transprob(U, H, energy, baseline)
end

"""
$(SIGNATURES)

Calculate the transistion probabilities between the neutrino flavours

# Arguments
- `osc_params::OscillationParameters`:  Oscillation parameters
- `energy`:                             Baseline [km]
- `baseline`:                           Energy [GeV]

"""
function transprob(osc_params::OscillationParameters, flavors::Pair{T, T}, energy, baseline) where {T <: Union{NeutrinoFlavor, Integer}}
    fromflavor = Int(first(flavors))
    toflavor = Int(last(flavors))
    transprob(osc, energy, baseline)[fromflavor, toflavor]
end

"""
$(SIGNATURES)

Returns the number of CP violating phases at given number of neutrino types

# Arguments
- `n`: number of neutrino types in the supposed model

# Examples
```julia-repl
julia> cpphases(3)
1
```
"""
function cpphases(n)
    n < 1 && return 0
    div((n - 1) * (n - 2), 2 )
end

"""
$(SIGNATURES)

Returns the number of mixing angles at given number of neutrino types

# Arguments
- `n`: number of neutrino types in the supposed model

# Examples
```julia-repl
julia> mixingangles(3)
3
```
"""
function mixingangles(n)
    n < 1 && return 0
    div(n * (n - 1), 2)
end
