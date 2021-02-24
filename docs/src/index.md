![](https://github.com/KM3NeT/Neurthino.jl/raw/master/docs/src/assets/neurthino.png)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KM3NeT.github.io/Neurthino.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KM3NeT.github.io/Neurthino.jl/dev)
[![Build Status](https://github.com/KM3NeT/Neurthino.jl/workflows/CI/badge.svg)](https://github.com/KM3NeT/Neurthino.jl/actions)
[![Codecov](https://codecov.io/gh/KM3NeT/Neurthino.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/KM3NeT/Neurthino.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3933364.svg)](https://doi.org/10.5281/zenodo.3933364)


# Neurthino.jl

**The package is currently under development**

**Neurthino.jl** is a package for calculating neutrino oscillation probabilities.
The main focus of the package lies on atmospheric neutrino flux and the neutrino
propagation through earth.

### Basic Usage
First of all the basic vacuum properties have to be defined by creating a
`OscillationParameters` struct with fixed number of neutrino flavours of the 
considered model:

```julia
julia> using Neurthino

julia> osc = OscillationParameters(3);
```

The values of the mixing angles (`setθ!`), mass squared differences (`setΔm²`)
and CP phases (`setδ!`) are initialised to 0 and have to be set individually:

```
julia> setθ!(osc, 1=>2, 0.59);

julia> setθ!(osc, 1=>3, 0.15);

julia> setθ!(osc, 2=>3, 0.84);

julia> setδ!(osc, 1=>3, 3.86);
```

The mass squared differences are defined as <img src="https://render.githubusercontent.com/render/math?math=\Delta_{ij}=m_i^2-m_j^2"> and
within the package the convention <img src="https://render.githubusercontent.com/render/math?math=\forall%20i%3Cj:m_i%3Cm_j"> is kept.

```
julia> setΔm²!(osc, 2=>3, -2.523e-3);

julia> setΔm²!(osc, 1=>2, -7.39e-5);
```

These oscillation parameters can now be used to calculate the oscillation
probabilities between the flavour states: 

```
julia> p = Pνν(osc, 1, 10000)
4-dimensional AxisArray{Float64,4,...} with axes:
    :Energy, [1.0]
    :Baseline, [10000.0]
    :InitFlav, NeutrinoFlavour[Electron, Muon, Tau]
    :FinalFlav, NeutrinoFlavour[Electron, Muon, Tau]
And data, a 1×1×3×3 Array{Float64,4}:
[:, :, 1, 1] =
 0.4027032498649196

[:, :, 2, 1] =
 0.10020685077639596

[:, :, 3, 1] =
 0.4970898993586843

[:, :, 1, 2] =
 0.24862049503031283

[:, :, 2, 2] =
 0.49124539964767405

[:, :, 3, 2] =
 0.2601341053220131

[:, :, 1, 3] =
 0.34867625510476735

[:, :, 2, 3] =
 0.4085477495759299

[:, :, 3, 3] =
 0.24277599531930263
```

The output is an `AxisArray` which provides intuitive indexing, e.g.
for P(νμ→ντ) at the given energy and baseline:

```
julia> p[Energy=1, Baseline=1, InitFlav=Muon, FinalFlav=Tau]
0.4085477495759299
```

The probabilities are calculated based on the transition matrix 
(the so-called PMNS-Matrix) between flavour and mass eigenstates,
as well as the Hamiltonian in the mass eigenbasis. In order to calculating these 
just once, the `Pνν` function can be utilised in the following way:

```
julia> U = PMNSMatrix(osc)
3×3 Array{Complex,2}:
   0.82161+0.0im         0.550114+0.0im        -0.112505+0.0983582im
 -0.301737+0.0608595im   0.601232+0.0407488im   0.736282+0.0im      
  0.476688+0.0545516im  -0.576975+0.0365253im   0.659968+0.0im

julia> H = Hamiltonian(osc)
3-element Array{Complex{Float64},1}:
 -0.0008902666666666667 + 0.0im
 -0.0008163666666666667 + 0.0im
  0.0017066333333333333 + 0.0im

julia> Pνν(U, H, 1, 10000)
4-dimensional AxisArray{Float64,4,...} with axes:
    :Energy, [1.0]
    :Baseline, [10000.0]
    :InitFlav, NeutrinoFlavour[Electron, Muon, Tau]
    :FinalFlav, NeutrinoFlavour[Electron, Muon, Tau]
And data, a 1×1×3×3 Array{Float64,4}:
[:, :, 1, 1] =
 0.4027032498649196

[:, :, 2, 1] =
 0.10020685077639596

[:, :, 3, 1] =
 0.4970898993586843

[:, :, 1, 2] =
 0.24862049503031283

[:, :, 2, 2] =
 0.49124539964767405

[:, :, 3, 2] =
 0.2601341053220131

[:, :, 1, 3] =
 0.34867625510476735

[:, :, 2, 3] =
 0.4085477495759299

[:, :, 3, 3] =
 0.24277599531930263
```

The neutrino propagation through matter can be calculated in two different ways.

For **homogeneous matter with a fixed density**, a modified PMNS-Matrix
and Hamiltonian can be determined and passed into `Pνν`, just like for
oscillations in vacuum. In order to determine the modified PMNS-Matrix and
Hamiltonian the neutrino energy and the matter density are required: 

```
julia> U_mat, H_mat = MatterOscillationMatrices(U, H, 1, 13);

julia> H_mat
3-element Array{Complex{Float64},1}:
 -0.0008318913110425209 + 1.402405683460369e-20im 
  0.0009755131119987117 + 4.535370547007611e-21im 
  0.0018408194174787428 - 2.1947559170628515e-20im

julia> U_mat
3×3 Array{Complex{Float64},2}:
  0.0358018-0.000158113im  0.970863+0.0im       -0.178275+0.156083im
 -0.662778+0.00661213im    0.157174+0.116074im   0.722845+0.0im
  0.74793+0.0im            0.0917808+0.104043im  0.649115-0.00104331im
```

The oscillation probabilities using the `Pνν` function, as described above:

```
julia> Pνν(U_mat, H_mat, 1, 10000)
4-dimensional AxisArray{Float64,4,...} with axes:
    :Energy, [1]
    :Baseline, [10000]
    :InitFlav, NeutrinoFlavour[Electron, Muon, Tau]
    :FinalFlav, NeutrinoFlavour[Electron, Muon, Tau]
And data, a 1×1×3×3 Array{Float64,4}:
[:, :, 1, 1] =
 0.8342726992356658

[:, :, 2, 1] =
 0.10810503748752147

[:, :, 3, 1] =
 0.057622263276811116

[:, :, 1, 2] =
 0.08277404351083649

[:, :, 2, 2] =
 0.052276251315704

[:, :, 3, 2] =
 0.8649497051734536

[:, :, 1, 3] =
 0.08295325725349613

[:, :, 2, 3] =
 0.839618711196769

[:, :, 3, 3] =
 0.07742803154973019

```

The second option is suitable for scenarios with more **complex paths** with
sections of different densities. An example is shown in the next chapter, where
we propagate neutrinos through the earth.

### Neutrino propagation through the Earth

The `Neurthino.jl` package also includes features for the neutrino oscillation probabilities
through the Earth, i.e. it contains functions for generating a neutrino path based on the
PREM model. In the following example a neutrino oscillogram with a resolution of 200x200 bins
is determined. The zenith angles for up going neutrinos (cos(θ)ϵ[-1,0]) and 
subsequently the neutrino paths are generated first:

```
julia> zenith = acos.(range(-1,stop=0,length=200));

julia> paths = Neurthino.prempath(zenith, 2.5, samples=100, discrete_densities=0:0.1:14);
```

The detector is assumed to be 2.5km under the earth's surface (a typical KM3NeT
detector block in the Mediterranean), which is a realistic scenario for
Water-Cherenkov-Detectors in sea or ice. Each path consists of 100 sections of
equal lengths while the matter density is taken from the PREM model.
If a vector of densities is passed as `discrete_densities`, the values are 
clipped to the closest value.

```
julia> energies = 10 .^ range(0, stop=2, length=200);

julia> prob = Pνν(U, H, energies, paths);
```
The returned array `prob` is again of type `AxisArray` with an axis `Path` for the path index (instead of the `Baseline` axis).
P(νe&#8594;νe) is determined by `prob[InitFlav=Electron, FinalFlav=Electron]`, which can be visualised by a `heatmap`:<br />
![](https://github.com/KM3NeT/Neurthino.jl/raw/master/docs/src/assets/earth_prob_elel.png) <br />
and for P(νμ&#8594;νμ) or `prob[InitFlav=Muon, FinalFlav=Muon]`:<br />
![](https://github.com/KM3NeT/Neurthino.jl/raw/master/docs/src/assets/earth_prob_mumu.png)
<!-- ```@index -->
<!-- ``` -->
<!--  -->
<!-- ```@autodocs -->
<!-- Modules = [Neurthino] -->
<!-- ``` -->
