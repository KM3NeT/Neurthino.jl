struct Path{T}
    density::Vector{T}
    baseline::Vector{T}
end

Path(density::Float64, baseline::Float64) = Path{Float64}([density],[baseline])

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
- `energy`: Neutrino energy [GeV]
- `density`: Matter density in g*cm^-3 
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino

"""
function MatterOscillationMatrices(osc_vacuum::OscillationParameters, energy, density; zoa=0.5, anti=false)
    H_vacuum = Diagonal(Hamiltonian(osc_vacuum)) 
    U_vacuum = PMNSMatrix(osc_vacuum)
    H_eff = convert(Array{ComplexF64}, U_vacuum * Diagonal{ComplexF64}(H_vacuum) * adjoint(U_vacuum))
    return MatterOscillationMatrices(H_eff, energy, density; zoa=zoa, anti=anti)
end

"""
$(SIGNATURES)

# Arguments
- `osc_vacuum::OscillationParameters`: Vacuum oscillation parameters
- `energy`: Neutrino energy [GeV]
- `paths`: Neutrino paths
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function transprob(osc_vacuum::OscillationParameters, energy, paths::Vector{Path{Float64}}; zoa=0.5, anti=false)
    # TODO: attach U_vac and H_vac to the oscillation parameters, so that it's
    # only calculated once and invalidated when any of the oscillation parameters
    # are changed
    U_vac = PMNSMatrix(osc_vacuum)
    H_vac = Hamiltonian(osc_vacuum)
    transprob(U_vac, H_vac, energy, densities, baselines; zoa=zoa, anti=anti)
end


"""
$(SIGNATURES)

# Arguments
- `U`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `energies`: Neutrino energies [GeV]
- `paths`: Neutrino paths
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function transprob(U, H, energies, paths::Vector{Path{Float64}}; zoa=0.5, anti=false)
    H_eff = U * Diagonal{ComplexF64}(H) * adjoint(U)
    A = zeros(ComplexF64, length(energies), length(paths), size(U)...)
    cache_size = length(energies) * sum(map(x->length(x.density), paths)) 
    lru = LRU{Tuple{Float64, Float64},
              Tuple{Array{ComplexF64,2}, Vector{ComplexF64}}}(maxsize=cache_size)
    for k in 1:length(energies)
        @inbounds E = energies[k]
        for (l, p) in enumerate(paths)
            tmp = Matrix{ComplexF64}(1I, size(U))
            for (m,b) in enumerate(p.baseline)
                @inbounds ρ = p.density[m]
                U_mat, H_mat = get!(lru, (E, ρ)) do
                    MatterOscillationMatrices(copy(H_eff), E, ρ; zoa=zoa, anti=anti)
                end  
                tmp *= Neurthino._transprobampl(U_mat, H_mat, E, b)
            end
            @inbounds A[k, l,  :, :] = tmp        
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    P
end

"""
$(SIGNATURES)

# Arguments
- `osc_vacuum::OscillationParameters`: Vacuum oscillation parameters
- `energy`: Neutrino energy [GeV]
- `paths`: Neutrino paths
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
Pνν(osc_vacuum, energies, paths::Vector{Path{Float64}}; zoa=0.5, anti=false) = transprob(osc_vacuum, energies, paths; zoa=zoa, anti=anti)

"""
$(SIGNATURES)

# Arguments
- `U`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `energies`: Neutrino energies [GeV]
- `paths`: Neutrino paths
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
Pνν(U, H, energies, paths::Vector{Path{Float64}}; zoa=0.5, anti=false) = transprob(U, H, energies, paths; zoa=zoa, anti=anti)
