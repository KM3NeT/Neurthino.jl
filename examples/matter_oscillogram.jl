#!/usr/bin/env julia

println("Loading libraries")
using Neurthino
using Unitful
using Plots
using ProgressMeter


function main()
    println("Initialising oscillation parameters")
    osc = OscillationParameters(3);
    mixingangle!(osc, 1 => 2, 0.58);
    mixingangle!(osc, 1 => 3, 0.15);
    mixingangle!(osc, 2 => 3, 0.737);
    cpphase!(osc, 1 => 3, 5.34);
    masssquareddiff!(osc, 2 => 3, -2.457e-3);
    masssquareddiff!(osc, 1 => 3, -2.457e-3-7.5e-5);
    masssquareddiff!(osc, 1 => 2, -7.5e-5);

    U = PMNSMatrix(osc)
    H = Hamiltonian(osc)

    coszenith = range(-1; stop=0, length=100)
    zenith = acos.(coszenith)
    energies = 10 .^(range(0;stop=2,length=100))
    
    data_points = zeros((3,3,length(zenith), length(energies)))


    @showprogress 1 "Calculating transition probabilities for different zenith angles" for (i,z) in enumerate(zenith)
        sec, dens = Neurthino.prempath(z, 2.5, samples=50)
        P = Neurthino.mattertransprob(osc, energies, dens, sec)
        for (k,l) in Iterators.product(fill(1:3, 2)...)
            data_points[k,l,i,:] = map(x->x[k,l], P)
        end
    end
    
    plot_arr = Array{Plots.Plot{Plots.GRBackend},1}()
    for (k,l) in Iterators.product(fill(1:3, 2)...)
        p = heatmap(data_points[k,l,:,:],  yticks=(1:20:100, coszenith[1:10:100]), xticks=(1:40:200, energies[1:10:50]), clim=(0,1))
        push!(plot_arr, p)
    end
    files = Vector{AbstractString}()
    @showprogress 1 "Creating plots" for (idx, p) in enumerate(plot_arr)
        f = "$idx.pdf"
        push!(files, f)
        savefig(p, f)
    end
end

main()







