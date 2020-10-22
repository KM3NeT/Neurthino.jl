EARTH_RADIUS = 6371.0u"km"

function tracklength(zenith::Quantity, zposition::Quantity)
"""
$(SIGNATURES)

Total path length through earth from detector position

# Arguments
- `zenith::Quantity`: Zenith angle of the path with respect to the detector frame
- `zposition::Quantity` Distance below the surface of the Earth (positive value) 
"""
    if isapprox(uconvert(u"°", zenith).val, 180)
        return 2*EARTH_RADIUS - zposition
    elseif isapprox(uconvert(u"°", zenith).val, 0)
        return zposition
    end
    
    theta = pi*u"rad" - zenith
    zprime = EARTH_RADIUS - zposition
    sintheta = sin(theta)
    costheta = cos(theta)
    halfcoord = sqrt(EARTH_RADIUS^2 - (zprime * sintheta)^2)
    chordsection = zprime * costheta
    if isapprox(uconvert(u"°", zenith).val, 90)
        return halfcoord
    else
        return halfcoord + chordsection
    end
end

function prempath(zenith::Quantity, zposition::Quantity; samples=100)
"""
$(SIGNATURES)
# Arguments
- `zenith::Quantity`: Zenith angle of the path with respect to the detector frame
- `zposition::Quantity` Distance below the surface of the Earth (positive value) 
- `samples` The number of steps with equal distance
"""
    trklen = tracklength(zenith, zposition)
    x = Array(range(0.0; stop=uconvert(u"km",trklen).val, length=samples))
    sections = (x[2:end] - x[1:end-1])
    total_pathlen = 0.5 * (x[2:end] + x[1:end-1])
    zprime = uconvert(u"km", EARTH_RADIUS - zposition).val
    radii = map(x -> sqrt(zprime^2 + x^2 + 2*zprime*x*cos(zenith)), total_pathlen)
    densities = PREM.(radii)
    sections*u"km", densities*u"g/cm^3"
end

function mattertransprob(osc::OscillationParameters ,energy ,densities, baselines)
"""
$(SIGNATURES)
# Arguments
- `osc_params::OscillationParameters`: Oscillation parameters
- `energy`: Neutrino energy
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
