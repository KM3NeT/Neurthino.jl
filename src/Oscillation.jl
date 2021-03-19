@enum NeutrinoFlavour begin
  Electron = 1
  Muon = 2
  Tau = 3
end


struct OscillationParameters{T}
    mixing_angles::SparseMatrixCSC{T,<:Integer}
    mass_squared_diff::SparseMatrixCSC{T,<:Integer}
    cp_phases::SparseMatrixCSC{T,<:Integer}
    dim::Int64
    OscillationParameters(dim::Int64) = begin
        new{ComplexF64}(
                spzeros(dim, dim),
                spzeros(dim, dim),
                spzeros(dim, dim),
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
- `value<:Real` The value which should be applied to the oscillation parameters

"""
function mixingangle!(osc::OscillationParameters, indices::Pair{T, T}, value::S) where {T <: Integer, S <: Real}
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

Set a mixing angle of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `args::Tuple{Pair{<:Integer, <:Integer}, <:Real}`: The indices of the mixing angle

"""
function mixingangle!(osc::OscillationParameters, (args::Tuple{Pair{T, T}, S})...) where {T <: Integer, S<: Real}
    for a in args
        mixingangle!(osc, first(a), last(a))
    end
end

const setθ! = mixingangle!

function _mass_matrix_fully_determined(osc::OscillationParameters)
    I, J, _ = findnz(osc.mass_squared_diff)
    set_elements = collect(zip(I, J))
    indices = Set([first.(set_elements)..., last.(set_elements)...])
    (length(indices) >= osc.dim) & (length(set_elements) >= (osc.dim - 1))  
end

function _mass_matrix_overdetermined(osc::OscillationParameters)
    I, J, _ = findnz(osc.mass_squared_diff)
    set_elements = collect(zip(I, J))
    indices = Set([first.(set_elements)..., last.(set_elements)...])
    length(set_elements) >= length(indices) 
end

function Base.isvalid(osc::OscillationParameters)
    return _mass_matrix_fully_determined(osc)
end

function _completed_mass_matrix(osc::OscillationParameters)
    tmp = Matrix(osc.mass_squared_diff)
    tmp = tmp - transpose(tmp)
    if _mass_matrix_fully_determined(osc)
        I, J, _ = findnz(osc.mass_squared_diff)
        given_idx = collect(zip(I, J))
        wanted_idx = filter(x->(x[1] < x[2]) & (x ∉ given_idx), collect(Iterators.product(1:osc.dim, 1:osc.dim)))
        graph = SimpleDiGraph(map(x->x!=0.0, tmp))
        for (from, to) in wanted_idx
            path = a_star(graph, from, to)
            for edge in path
                tmp[from, to] += tmp[edge.src, edge.dst]
            end
        end
    else
        error("Mass squared differences not fully determined!")
    end
    UpperTriangular(tmp)
end

"""
$(SIGNATURES)

Set a mass squared difference of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Pair{<:Integer, <:Integer}`: The indices of the mass squared difference
- `value` The value which should be applied to the oscillation parameters

"""
function masssquareddiff!(osc::OscillationParameters, indices::Pair{T, T}, value::S) where {T <: Integer, S <: Number}
    fromidx = first(indices)
    toidx = last(indices)
    if fromidx < toidx
        osc.mass_squared_diff[fromidx, toidx] = value
    elseif fromidx == toidx
        error("Mass squared difference with equal index cannot be modified.")
    else
        osc.mass_squared_diff[toidx, fromidx] = -value
    end
    if _mass_matrix_overdetermined(osc)
        @warn "Mass squared difference fields (partially) overdetermined!"
    end
end

"""
$(SIGNATURES)

Set a mass squared difference of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `args::Tuple{Pair{<:Integer, <:Integer}, <:Number}`: Indices and values of the mass squared difference

"""
function masssquareddiff!(osc::OscillationParameters, (args::Tuple{Pair{<:Integer, <:Integer}, <:Number})...)
    for a in args
        masssquareddiff!(osc, first(a), last(a))
    end
end

const setΔm²! = masssquareddiff!

"""
$(SIGNATURES)

Set a CP phase of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `indices::Pair{<:Integer, <:Integer}`: The indices of the mass difference
- `value` The value which should be applied to the oscillation parameters

"""
function cpphase!(osc::OscillationParameters, indices::Pair{T, T}, value::S) where {T <: Integer, S <: Real}
    fromidx = first(indices)
    toidx = last(indices)
    if fromidx < toidx
        osc.cp_phases[fromidx, toidx] = value
    else
        osc.cp_phases[toidx, fromidx] = value
    end
end

"""
$(SIGNATURES)

Set a CP phase of an oscillation parameters struct

# Arguments
- `osc::OscillationParameters`: Oscillation parameters 
- `args::Tuple{Pair{<:Integer, <:Integer}, <:Number}`: Indices and values of the CP phase

"""
function cpphase!(osc::OscillationParameters, (args::Tuple{Pair{T, T}, S})...) where {T <: Integer, S <: Real}
    for a in args
        cpphase!(osc, first(a), last(a))
    end
end

const setδ! = cpphase!

function PMNSMatrix(osc_params::OscillationParameters; anti=false)
"""
$(SIGNATURES)

Create rotation matrix (PMNS) based on the given oscillation parameters

# Arguments
- `osc_params::OscillationParameters`: Oscillation parameters
- `anti`: Is anti neutrino

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
            if anti
                cp_term = conj(cp_term)
            end
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
    full_mass_squared_matrix = _completed_mass_matrix(osc_params)
    H = zeros(ComplexF64, osc_params.dim)
    for i in 1:osc_params.dim
        for j in 1:osc_params.dim
            if i < j
                H[i] += full_mass_squared_matrix[i,j]
            elseif j < i
                H[i] -= full_mass_squared_matrix[j,i]
            end
        end
        H[i] += 1im * lambda[i]
    end
    H /= osc_params.dim
    H
end


function _oscprobampl(U, H, energy, baseline)  
    H_diag = 2.5338653580781976 * Diagonal{ComplexF64}(H) * baseline / energy 
    U * exp(1im * H_diag) * adjoint(U)
end

function _make_flavour_range(size::Integer)
    if size <= 3
        return NeutrinoFlavour.(1:size)
    else
        return [NeutrinoFlavour.(1:3)..., 4:size...]
    end
end


"""
$(SIGNATURES)

Calculate the transistion probabilities between the neutrino flavours

# Arguments
- `U`:          PMNS Matrix
- `H`:          Hamiltonian
- `energy`:     Energies [GeV]
- `baseline`:   Baselines [km]

"""
function oscprob(U, H, energy::Vector{T}, baseline::Vector{S}) where {T,S <: Real}
    s = (size(U)..., length(energy), length(baseline))
    combinations = collect(Iterators.product(energy, baseline))
    tmp = map(x->abs.(_oscprobampl(U, H, first(x), last(x))).^2, combinations)
    P = reshape(hcat(collect(Iterators.flatten(tmp))), s...)
    P = permutedims(P, (3,4,1,2))
    flavrange = _make_flavour_range(first(size(U)))
    AxisArray(P; Energy=energy, Baseline=baseline, InitFlav=flavrange, FinalFlav=flavrange)
end

const oscprob(U, H, energy::T, baseline::Vector{S}) where {S,T <: Real} = oscprob(U, H, [energy], baseline)
const oscprob(U, H, energy, baseline::T) where {T <: Real} = oscprob(U, H, energy, [baseline])

"""
$(SIGNATURES)

Calculate the transistion probabilities between the neutrino flavours

# Arguments
- `osc_params::OscillationParameters`:  Oscillation parameters
- `energy`:                             Energy [GeV]
- `baseline`:                           Baseline [km]
- `anti`:                               Is anti neutrino

"""
function oscprob(osc_params::OscillationParameters, energy, baseline; anti=false)  
    H = Hamiltonian(osc_params)
    U = PMNSMatrix(osc_params; anti=anti)
    Pνν(U, H, energy, baseline)
end

const Pνν = oscprob

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
