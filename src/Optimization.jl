function interpolatematteroscillationmatrices(osc::OscillationParameters, densities, energies)
    H = Hamiltonian(osc)
    U = PMNSMatrix(osc)
    interpolatematteroscillationmatrices(U, H, densities, energies)
end

function interpolatematteroscillationmatrices(U, H, densities, energies)
    knots = collect(Iterators.product(densities, energies))
    matrices = map(x->MatterOscillationMatrices(U, H, x...), knots)
    H_shape = (size(H)..., size(densities)..., size(energies)...)
    U_shape = (size(U)..., size(densities)..., size(energies)...)
    H_samples = reshape(collect(Iterators.flatten(map(x->abs.(x[2]), matrices))), H_shape)
    U_samples = reshape(collect(Iterators.flatten(map(x->x[1], matrices))), U_shape)
    function relative_phases(A)
        phases = angle.(A)
        retval = zeros(Float64, size(A)...)
        for i in 1:size(A)[1]
            retval[:,i] = phases[:,i] .- phases[i,i]
        end
        retval
    end
    U_phases = mapslices(relative_phases, U_samples, dims=[1,2])
    abs_itp(d,e) = mapslices(x->CubicSplineInterpolation((densities, energies), abs.(x))(d,e), U_samples, dims=[3,4])
    phase_itp(d,e) = mapslices(x->CubicSplineInterpolation((densities, energies), x)(d,e), U_phases, dims=[3,4])
    U_itp(d,e) = (abs_itp(d,e) .* exp.(-1im .* phase_itp(d,e)))[:,:,1,1]
    H_itp(d,e) = mapslices(x->CubicSplineInterpolation((densities, energies), abs.(x))(d,e), H_samples, dims=[2,3])[:,1,1]
    U_itp, H_itp
end
