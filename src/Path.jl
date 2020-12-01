function mattertransprob(osc::OscillationParameters, energy, densities, baselines)
"""
$(SIGNATURES)

# Arguments
- `osc_params::OscillationParameters`: Vacuum oscillation parameters
- `energy`: Neutrino energy [GeV]
- `densities`: Matter densities of the traversed path [g/cm^3]
- `baselines`: Path section lengths [km]
"""
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
            H_mat, U_mat = MatterOscillationMatrices(U, H, densities[i], energy[n])
            H_tmp = 2.534 * Diagonal(H_mat) * b / energy[n]
            A[n] *= U_mat * exp(-1im * H_tmp) * adjoint(U_mat)
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    P
end
