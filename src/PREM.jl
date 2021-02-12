const EARTH_RADIUS = 6371.0  # km

struct Layer
    name::AbstractString
    r_min
    r_max
    density::Polynomial
end


struct EarthModel
    layers::Vector{Layer}
end


PREM = EarthModel(
    [
        Layer("Inner Core", 0.0, 1221.5, Polynomial([13.0885, 0.0, -8.8381])),
        Layer("Outer Core", 1221.5, 3480.0, Polynomial([12.5815, -1.2638, -3.6426, -5.5281])),
        Layer("Lower Mantle", 3480.0, 5701.0, Polynomial([7.9565, -6.4761, 5.5283, -3.0807])),
        Layer("Transition Zone 1", 5701.0, 5771.0, Polynomial([5.3197, -1.4836])),
        Layer("Transition Zone 2", 5771.0, 5971.0, Polynomial([11.2494, -8.0298])),
        Layer("Transition Zone 3", 5971.0, 6151.0, Polynomial([7.1089, -3.8045])),
        Layer("LVZ", 6151.0, 6291.0, Polynomial([2.6910, 0.6924])),
        Layer("LID", 6291.0, 6346.6, Polynomial([2.6910, 0.6924])),
        Layer("Crust 1", 6346.6, 6356.0, Polynomial([2.9])),
        Layer("Crust 2", 6356.0, 6368.0, Polynomial([2.6])),
        Layer("Ocean", 6368.0, 6371.0, Polynomial([1.020]))
    ]
)


function (m::EarthModel)(r)
    for layer ∈ m.layers
        if r < layer.r_max && r >= layer.r_min
            return layer.density(r / m.layers[end].r_max)
        end
    end
end


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

function prempath(zenith::Vector{Float64}, zposition; samples=100, discrete_densities=nothing)
"""
$(SIGNATURES)

# Arguments
- `zenith::Vector{Float64}`: Zenith angles of the paths with respect to the detector frame [rad]
- `zposition::Float64` Distance below the surface of the Earth (positive value) [km]
- `samples` The number of steps with equal distance
- `discrete_densities` List of density values to be used for discretization
"""
    map(z->prempath(z, zposition; samples=samples, discrete_densities=discrete_densities), zenith)
end

function prempath(zenith::T, zposition; samples=100, discrete_densities=nothing) where {T <: Number}
"""
$(SIGNATURES)

# Arguments
- `zenith::Float64`: Zenith angle of the path with respect to the detector frame [rad]
- `zposition::Float64` Distance below the surface of the Earth (positive value) [km]
- `samples` The number of steps with equal distance
"""
    trklen = tracklength(zenith, zposition)
    x = Array(range(0.0; stop=trklen, length=samples))
    sections = (x[2:end] - x[1:end-1])
    total_pathlen = 0.5 * (x[2:end] + x[1:end-1])
    zprime = EARTH_RADIUS - zposition
    radii = map(x -> sqrt(zprime^2 + x^2 + 2*zprime*x*cos(zenith)), total_pathlen)
    densities = PREM.(radii)
    if !isnothing(discrete_densities)
        idx = map(d->searchsortedfirst(discrete_densities,d), densities)
        densities = map(i->discrete_densities[i], idx)
    end
    Path(densities, sections)
end
