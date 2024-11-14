#!/usr/bin/env julia

using Neurthino
using Unitful
using Plots
using LaTeXStrings

function example_matter_oscillogram()
    # run the example_matter_oscillogram

    # define oscillation parameters with 3 neutrino flavours
    osc = OscillationParameters(3);

    # define values of the mixing angles, cpp phase and mass square differences
    mixingangle!(osc, 1 => 2, 0.58);
    mixingangle!(osc, 1 => 3, 0.15);
    mixingangle!(osc, 2 => 3, 0.737);
    cpphase!(osc, 1 => 3, 5.34);
    masssquareddiff!(osc, 2 => 3, -2.457e-3);
    masssquareddiff!(osc, 1 => 2, -7.5e-5);

    # compute the PMNS matrix and Hamiltonian before starting the loop
    U = PMNSMatrix(osc)
    H = Hamiltonian(osc)

    # the zenith and energy domain to compute
    coszenith = range(-1; stop=0, length=100)
    zenith = acos.(coszenith)
    energies = 10 .^(range(0;stop=2,length=100))

    # initialise data array
    data_points = zeros((3,3,length(zenith), length(energies)))

    # loop over zenith directions
    for (i,z) in enumerate(zenith)

        # create a path at with zenith z and depth 2.5 km
        path = Neurthino.prempath(z, 2.5, samples=50)

        # compute the probabilities
        prob = Neurthino.oscprob(osc, energies, path)

        # save the results in data array
        for (k,l) in Iterators.product(fill(1:3, 2)...)
            data_points[k,l,i,:] = prob[:,1,k,l]
        end
    end

    # loop over the start and end flavours
    for (k, startFlav) in enumerate([L"\nu_e", L"\nu_\mu", L"\nu_\tau"])
        for (l, endFlav) in enumerate([L"\nu_e", L"\nu_\mu", L"\nu_\tau"])

            # create an oscillogram for every start- and endFlav combination
            p = heatmap(energies*u"GeV", coszenith, data_points[k,l,:,:], xscale=:log10,
                xlabel="Energy", ylabel=L"\cos(\theta_z)", title=latexstring(startFlav, L"\rightarrow", endFlav),
                clim=(0,1));

            # save each plot
            # f = "nu_$(k)_nu_$(l).pdf"
            # savefig(p, f)
            # which is skipped for testing purposes
        end
    end
end