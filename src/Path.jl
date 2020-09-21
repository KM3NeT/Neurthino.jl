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
