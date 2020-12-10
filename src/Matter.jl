"""
$(SIGNATURES)

Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `U`: Vacuum PMNSMatrix
- `H`: Vacuum Hamiltonian
- `density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(U, H, energy, density; zoa=0.5, anti=false)
    H_eff = U * Diagonal{Complex}(H) * adjoint(U)
    H_eff = H_eff
    A = sqrt(2) * G_F * N_A * zoa * density
    if anti
        H_eff[1,1] -= A * (2 * energy * 1e9)
    else
        H_eff[1,1] += A * (2 * energy * 1e9)
    end
    U_matter = eigvecs(H_eff)
    H_matter = eigvals(H_eff)
    return U_matter, H_matter
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
    return MatterOscillationMatrices(U_vacuum, H_vacuum, energy, density)
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
    U_vac = PMNSMatrix(osc)
    H_vac = Hamiltonian(osc)
    mattertransprob(U_vac, H_vac, energy, densities, baselines)
end


function mattertransprob(U, H, energy, densities, baselines)
"""
$(SIGNATURES)

# Arguments
- `U`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `energy`: Neutrino energy [GeV]
- `densities`: Matter densities of the traversed path [g/cm^3]
- `baselines`: Path section lengths [km]
"""
    A = fill(Matrix{Complex}(1I, size(U)), length(energy))
    for n in 1:length(energy)
        for (i,b) in enumerate(baselines)
            U_mat, H_mat = MatterOscillationMatrices(U, H, energy[n], densities[i])
            A[n] *= Neurthino._transprobampl(U_mat, H_mat, energy[n], b)
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    P
end
