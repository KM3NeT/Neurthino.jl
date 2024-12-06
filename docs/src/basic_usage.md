# Basic usage

## Setting up the Oscillation Parameters
We start by creating an `OscillationParameters` struct in vacuum with 3 neurtrino flavours:
```@example 1
using Neurthino

osc = OscillationParameters(3);
nothing # hide
```

The mixing angles, CP-phase and mass squared differences are all initialised to be zero. We can set them individually by:
```@example 1
mixingangle!(osc, 1 => 2, 0.58);
mixingangle!(osc, 1 => 3, 0.15);
mixingangle!(osc, 2 => 3, 0.737);

cpphase!(osc, 1 => 3, 5.34);

masssquareddiff!(osc, 2 => 3, -2.457e-3);
masssquareddiff!(osc, 1 => 2, -7.5e-5);
nothing # hide
```
One could also use `setθ!`, `setδ!` and `setΔm²!`.

We can then use `oscprob`, or similar `Pνν`, to compute the oscillation probabilities for all 9 flavour changes. For an single energy of 1.5 GeV and baseline of 10000 km (we could also pass in a vector), we get:
```@example 1
prob = oscprob(osc, 1.5, 10000)
```
The output is an `AxisArray` which provides intuitive indexing, e.g.
for `prob(νμ→ντ)` at the given energy and baseline:
```@example 1
prob[Energy=1, Baseline=1, InitFlav=Muon, FinalFlav=Tau]
```
The probabilities are calculated based on the transition matrix
(the so-called PMNS-Matrix) between flavour and mass eigenstates,
as well as the Hamiltonian in the mass eigenbasis. We can compute PMNS matrix:
```@example 1
U = PMNSMatrix(osc)
```
and the Hamiltonian:
```@example 1
H = Hamiltonian(osc)
```
once, and avoid repeated calculations during a loop. Then, we can use `oscprob(U, H, 1, 10000)` to get the same results as above.

!!! note "Oscillations in vacuum"
    Up to this point, the `osc`, `U` and `H` are all defined in vacuum. Using multiple dispatch, the same function `oscprob` will be used to get the oscillation probabilities along a path in matter.

## Neutrino oscillations in water
`Neurthino.jl` comes with the PREM model, a layered Earth model, which allows us to create a neutrino path based on just the depth and zenith angle. Instead of moving through a vacuum, we want to model a Water-Cherenkov detector in sea or ice. For this example, we place the detector at 2km depth:
```@example 1
waterdepth_km = 2;
nothing # hide
```
The zenith angles for up going neutrinos  and
subsequently the neutrino paths are generated first:

We first set the desired zenith angles (cos(θ)ϵ[-1,0]) for up going neutrinos and with an energy range of [0, 2] GeV:
```@example 1
coszenith = range(-1; stop=0, length=100);
zenith = acos.(coszenith);

energies = 10 .^(range(0;stop=2,length=100));
nothing # hide
```


With the zenith angles defined, we will create a loop that computes the oscillation probabilities for each angle.

First we pre-allocate the data array
```@example 1
data_points = zeros((3,3,length(zenith), length(energies)));
nothing # hide
```

We can pass the the matrices `U` and `H` directly to avoid repeating their computation. Although `U` and `H` were defined in vacuum, the call to `oscprob` will automatically convert them using the densities defined by `path`. Let's start the loop over zeniths:
```@example 1
for (i,z) in enumerate(zenith)

    # create a path with zenith z and waterdepth_km
    path = Neurthino.prempath(z, waterdepth_km, samples=50);

    # compute the probabilities for all flavours
    prob = oscprob(U, H, energies, path);

    # save the results in data array
    for (k,l) in Iterators.product(fill(1:3, 2)...)
        data_points[k,l,i,:] = prob[:,1,k,l];
    end
end
```
Here, we used `Neurthino.prempath` to create a `path`, which defines the depths and the densities according to the PREM Earth model. The layered model is defined by `Neurthino.PREM`. The first 3km is water.

Then we call `oscprob` to compute the oscillation probabilities of all neutrinos on that path. Note that we can use both `oscprob(osc, energies, path)` and `oscprob(U, H, energies, path)` (latter is prefered during loops to avoid recalculating `U` and `H` in each iteration). Lastly, we save the probabilities into our pre-allocated array.

Now we have our data, let's visualize the oscillogram! First we import the packages we will use:
```@example 1
using Unitful
using Plots
using LaTeXStrings
```
We again loop over all flavours and create an heatmap for each combination:
```@example 1
for (k, startFlav) in enumerate([L"\nu_e", L"\nu_\mu", L"\nu_\tau"])
    for (l, endFlav) in enumerate([L"\nu_e", L"\nu_\mu", L"\nu_\tau"])

        # create an oscillogram for every start- and endFlav combination
        p = heatmap(energies*u"GeV", coszenith, data_points[k,l,:,:], xscale=:log10,
            xlabel="Energy", ylabel=L"\cos(\theta_z)", title=latexstring(startFlav, L"\rightarrow", endFlav),
            clim=(0,1));

        # save each plot
        filename = "water_nu_$(k)_nu_$(l).svg"
        savefig(p, filename)

    end
    # instead of showing the output for each flavour, we break here
    break
end
```

For the electron neutrino oscillations, we then get:
![water_oscillogram](water_nu_1_nu_1.svg)

![water_oscillogram](water_nu_1_nu_2.svg)

![water_oscillogram](water_nu_1_nu_3.svg)


## Vectorize operations
The for-loops can be avoided by vectorizing the calculations. We first define all the paths based on all zenith angles and the waterdepth:
```@example 1
paths = Neurthino.prempath(zenith, waterdepth_km, samples=100);
nothing # hide
```

Then we compute all probabilites at once:
```@example 1
prob = oscprob(U, H, energies, paths)
```
The returned array `prob` is again of type `AxisArray` with an axis `:Energy` and now `:Path` for the path index (instead of the `:Baseline` axis).