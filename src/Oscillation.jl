struct OscillationParameters{T, N}
    mixing_angles::UnitUpperTriangular{T, <:MMatrix{N, N, T}}
    mass_squared_diff::UnitUpperTriangular{T, <:MMatrix{N, N, T}}
    cp_phases::UnitUpperTriangular{T, <:MMatrix{N, N, T}}

    OscillationParameters(dim::Integer) = begin
        new{Float64, dim}(
                UnitUpperTriangular(@MMatrix zeros(dim, dim)),
                UnitUpperTriangular(@MMatrix zeros(dim, dim)),
                UnitUpperTriangular(@MMatrix zeros(dim, dim)))
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
- `indices::Tuple{<:Integer, <:Integer}`: The indices of the mixing angle
- `value` The value which should be applied to the oscillation parameters

"""
function mixingangle!(osc::OscillationParameters, indices::Tuple{T, T}, value) where {T <: Integer}
    if indices[1] < indices[2]
        osc.mixing_angles[indices[1], indices[2]] = value
    else
        osc.mixing_angles[indices[2], indices[1]] = value
    end
end

"""
$(SIGNATURES)

Set a mass squared difference of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Tuple{<:Integer, <:Integer}`: The indices of the mass squared difference
- `value` The value which should be applied to the oscillation parameters

"""
function masssquareddiff!(osc::OscillationParameters, indices::Tuple{T, T}, value) where {T <: Integer}
    if indices[1] < indices[2]
        osc.mass_squared_diff[indices[1], indices[2]] = value
    else
        osc.mass_squared_diff[indices[2], indices[1]] = -value
    end
end

"""
$(SIGNATURES)

Set a CP phase of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Tuple{<:Integer, <:Integer}`: The indices of the mass difference
- `value` The value which should be applied to the oscillation parameters

"""
function cpphase!(osc::OscillationParameters, indices::Tuple{T, T}, value) where {T <: Integer}
    if indices[1] < indices[2]
        osc.cp_phases[indices[1], indices[2]] = value
    else
        osc.cp_phases[indices[2], indices[1]] = value
    end
end

"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `osc_vacuum::OscillationParameters`: Oscillation parameters in vacuum
- `matter_density`: Matter density in g*cm^-3 

"""
function MatterOscillationMatrices(osc_vacuum::OscillationParameters, matter_density)
    H_vacuum = Diagonal(Hamiltonian(osc_vacuum)) 
    U_vacuum = PMNSMatrix(osc_vacuum)
    return MatterOscillationMatrices(U_vacuum, H_vacuum, matter_density)
end

"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `U`: Vacuum PMNSMatrix
- `H`: Vacuum Hamiltonian
- `matter_density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(U, H, matter_density, energy; zoa=0.5, anti=false)
    H_eff = U * Diagonal{Complex}(H) * adjoint(U)
    H_eff = H_eff
    A = sqrt(2) * G_F * N_A * zoa * matter_density
    if anti
        H_eff[1,1] -= A * (2 * energy * 1e9)
    else
        H_eff[1,1] += A * (2 * energy * 1e9)
    end
    U_matter = eigvecs(H_eff)
    H_matter = eigvals(H_eff)
    return H_matter, U_matter
end

function PMNSMatrix(osc_params::OscillationParameters)
"""
$(SIGNATURES)

Create rotation matrix (PMNS) based on the given oscillation parameters

# Arguments
- `osc_params::OscillationParameters`: Oscillation parameters

"""
    dim = size(osc_params.mixing_angles)[1]
    pmns = Matrix{Complex}(1.0I, dim, dim) 
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
    dim = size(osc_params.mixing_angles)[1]
    Hamiltonian(osc_params, zeros(Float64, dim))
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
    dim = size(osc_params.mixing_angles)[1]
    H = zeros(Complex, dim)
    for i in 1:dim
        for j in 1:dim
            if i < j
                H[i] += osc_params.mass_squared_diff[i,j]
            elseif j < i
                H[i] -= osc_params.mass_squared_diff[j,i]
            end
        end
        H[i] += 1im * lambda[i]
    end
    H /= dim
    H
end


"""
$(SIGNATURES)

Calculate the transistion probability between the neutrino flavours

# Arguments
- `U`:          Unitary transition matrix
- `H`:          Energy eigenvalues
- `energy`:     Baseline [km]
- `baseline`:   Energy [GeV]

"""
function transprob(U, H, energy, baseline)  
    H_diag = 2.534 * baseline .* ( energy .+ Diagonal(H) ./ ( 2 * energy ) )
    A = U * exp(-1im * H_diag) * adjoint(U)
    P = abs.(A) .^ 2
end

"""
$(SIGNATURES)

Calculate the transistion probability between the neutrino flavours

# Arguments
- `osc_params::OscillationParameters`:  Oscillation parameters
- `energy`:                             Baseline [km]
- `baseline`:                           Energy [GeV]

"""
function transprob(osc_params::OscillationParameters, energy, baseline)  
    H = Hamiltonian(osc_params)
    U = PMNSMatrix(osc_params)
    H_diag = 2.534 * Diagonal(H) * baseline / energy 
    A = U * exp(-1im * H_diag) * adjoint(U)
    P = abs.(A) .^ 2
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
