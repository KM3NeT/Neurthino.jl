EARTH_RADIUS = 6371.0u"km"

function tracklength(zenith::Quantity, zposition::Quantity)
"""
$(SIGNATURES)

Total path lengths through earth from detector position

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
    x = Array(range(0, uconvert(u"km",trklen).val, length=samples))
    sections = (x[2:end] - x[1:end-1])
    total_pathlen = 0.5 * (x[2:end] + x[1:end-1])
    zprime = uconvert(u"km", EARTH_RADIUS - zposition).val
    radii = map(x -> sqrt(zprime^2 + x^2 + 2*zprime*x*cos(zenith)), total_pathlen)
    densities = PREM.(radii)
    sections*u"km", densities*u"g/cm^3"
end



