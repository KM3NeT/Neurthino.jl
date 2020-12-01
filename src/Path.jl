EARTH_RADIUS = 6371.0

function tracklength(zenith, zposition)
"""
$(SIGNATURES)

Total path length through earth from detector position

# Arguments
- `zenith`: Zenith angle of the path with respect to the detector frame [rad]
- `zposition` Distance below the surface of the Earth (positive value) [km]
"""
    if isapprox(zenith, pi)
        return 2*EARTH_RADIUS - zposition
    elseif isapprox(zenith, 0)
        return zposition
    end
    
    theta = pi - zenith
    zprime = EARTH_RADIUS - zposition
    sintheta = sin(theta)
    costheta = cos(theta)
    halfcoord = sqrt(EARTH_RADIUS^2 - (zprime * sintheta)^2)
    chordsection = zprime * costheta
    if isapprox(zenith, pi/2)
        return halfcoord
    else
        return halfcoord + chordsection
    end
end

function prempath(zenith, zposition; samples=100)
"""
$(SIGNATURES)

# Arguments
- `zenith::Quantity`: Zenith angle of the path with respect to the detector frame [rad]
- `zposition::Quantity` Distance below the surface of the Earth (positive value) [km]
- `samples` The number of steps with equal distance
"""
    trklen = tracklength(zenith, zposition)
    x = Array(range(0.0; stop=trklen, length=samples))
    sections = (x[2:end] - x[1:end-1])
    total_pathlen = 0.5 * (x[2:end] + x[1:end-1])
    zprime = EARTH_RADIUS - zposition
    radii = map(x -> sqrt(zprime^2 + x^2 + 2*zprime*x*cos(zenith)), total_pathlen)
    densities = PREM.(radii)
    sections, densities
end

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
    A = fill(Matrix{Complex}(1I, size(U_vac)), length(energy))
    for n in 1:length(energy)
        for (i,b) in enumerate(baselines)
            H_mat, U_mat = MatterOscillationMatrices(U_vac, H_vac, densities[i], energy[n])
            H_tmp = 2.534 * Diagonal(H_mat) * b / energy[n]
            A[n] *= U_mat * exp(-1im * H_tmp) * adjoint(U_mat)
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    P
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
