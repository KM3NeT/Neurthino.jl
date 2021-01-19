struct Path{T}
    density::Vector{T}
    baseline::Vector{T}
end

Base.iterate(p::Path, state=1) = state > length(p.density) ? nothing : ( (p.density[state], p.baseline[state]),  state+1)

Base.length(p::Path) = length(p.density)

"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `P`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(U, H, energy, density; zoa=0.5, anti=false)
    H_eff = convert(Array{ComplexF64}, U * Diagonal{Complex}(H) * adjoint(U))
    MatterOscillationMatrices(H_eff, energy, density; zoa=zoa, anti=anti)
end


"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `H_eff`: Effective Matter Hamiltonian
- `density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(H_eff, energy, density; zoa=0.5, anti=false)
    A = sqrt(2) * G_F * N_A * zoa * density
    a = H_eff[1,1]
    if anti
        H_eff[1,1] -= A * (2 * energy * 1e9)
    else
        H_eff[1,1] += A * (2 * energy * 1e9)
    end
    tmp = eigen(H_eff)
    return tmp.vectors, tmp.values
end

"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `osc_vacuum::OscillationParameters`: Oscillation parameters in vacuum
- `density`: Matter density in g*cm^-3 

"""
function MatterOscillationMatrices(osc_vacuum::OscillationParameters, energy, density)
    H_vacuum = Diagonal(Hamiltonian(osc_vacuum)) 
    U_vacuum = PMNSMatrix(osc_vacuum)
    H_eff = convert(Array{ComplexF64}, U_vacuum * Diagonal{Complex}(H_vacuum) * adjoint(U_vacuum))
    return MatterOscillationMatrices(H_eff, energy, density)
end

"""
$(SIGNATURES)

# Arguments
- `osc_params::OscillationParameters`: Vacuum oscillation parameters
- `energy`: Neutrino energy [GeV]
- `densities`: Matter densities of the traversed path [g/cm^3]
- `baselines`: Path section lengths [km]
"""
function mattertransprob(osc::OscillationParameters, energy, densities, baselines)
    # TODO: attach U_vac and H_vac to the oscillation parameters, so that it's
    # only calculated once and invalidated when any of the oscillation parameters
    # are changed
    U_vac = PMNSMatrix(osc)
    H_vac = Hamiltonian(osc)
    mattertransprob(U_vac, H_vac, energy, densities, baselines)
end


"""
$(SIGNATURES)

# Arguments
- `U`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `energies`: Neutrino energies [GeV]
- `paths`: Neutrino paths
"""
function mattertransprob(U, H, energies, paths::Array{Path{Float64},1})
    H_eff = U * Diagonal{Complex}(H) * adjoint(U)
    A = zeros(Complex, length(energies), length(paths), size(U)...)
    cache_size = length(energies) * sum(map(x->length(x.density), paths)) 
    lru = LRU{Tuple{Float64, Float64},
              Tuple{Array{Complex{Float64},2}, Vector{Complex{Float64}}}}(maxsize=cache_size)
    for k in 1:length(energies)
        @inbounds E = energies[k]
        for (l, p) in enumerate(paths)
            @inbounds A[k, l, :, :] = Matrix{ComplexF64}(1I, size(U))
            for (m,b) in enumerate(p.baseline)
                @inbounds ρ = p.density[m]
                U_mat, H_mat = get!(lru, (E, ρ)) do
                    MatterOscillationMatrices(copy(H_eff), E, ρ)
                end
                @inbounds A[k, l,  :, :] *= Neurthino._transprobampl(U_mat, H_mat, E, b)
            end
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    P
end
