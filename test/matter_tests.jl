using Neurthino

osc = OscillationParameters(3)
setθ!(osc, 1=>2, 0.59);
setθ!(osc, 1=>3, 0.15);
setθ!(osc, 2=>3, 0.84);
setδ!(osc, 1=>3, 3.86);
setΔm²!(osc, 2=>3, -2.523e-3);
setΔm²!(osc, 1=>2, -7.39e-5);
setδ!(osc, 1=>3, 234 * π / 180);

U_mat, H_mat = MatterOscillationMatrices(osc, 1, 13);
test_values = Pνν(U_mat, H_mat, 1, 10000)[Energy=1, Baseline=1]
@test test_values[1,1] ≈ 0.834 atol=0.01
@test test_values[1,2] ≈ 0.110 atol=0.01
@test test_values[1,3] ≈ 0.055 atol=0.01
@test test_values[2,1] ≈ 0.079 atol=0.01
@test test_values[2,2] ≈ 0.054 atol=0.01
@test test_values[2,3] ≈ 0.867 atol=0.01
@test test_values[3,1] ≈ 0.083 atol=0.01
@test test_values[3,2] ≈ 0.836 atol=0.01
@test test_values[3,3] ≈ 0.077 atol=0.01

p = Neurthino.Path([13],[10000])
E = [1,2]
test_values = Pνν(osc, E, p)
@test test_values[1,1,1,1] ≈ 0.834 atol=0.01
@test test_values[1,1,1,2] ≈ 0.110 atol=0.01
@test test_values[1,1,1,3] ≈ 0.055 atol=0.01
@test test_values[1,1,2,1] ≈ 0.079 atol=0.01
@test test_values[1,1,2,2] ≈ 0.054 atol=0.01
@test test_values[1,1,2,3] ≈ 0.867 atol=0.01
@test test_values[1,1,3,1] ≈ 0.083 atol=0.01
@test test_values[1,1,3,2] ≈ 0.836 atol=0.01
@test test_values[1,1,3,3] ≈ 0.077 atol=0.01


p = Neurthino.Path([13,13,13,13],[2500,2500,2500,2500])
E = 1

test_values = Pνν(osc, E, p)
@test test_values[1,1,1,1] ≈ 0.834 atol=0.01
@test test_values[1,1,1,2] ≈ 0.110 atol=0.01
@test test_values[1,1,1,3] ≈ 0.055 atol=0.01
@test test_values[1,1,2,1] ≈ 0.079 atol=0.01
@test test_values[1,1,2,2] ≈ 0.054 atol=0.01
@test test_values[1,1,2,3] ≈ 0.867 atol=0.01
@test test_values[1,1,3,1] ≈ 0.083 atol=0.01
@test test_values[1,1,3,2] ≈ 0.836 atol=0.01
@test test_values[1,1,3,3] ≈ 0.077 atol=0.01

@test test_values[Energy=1,Path=1,InitFlav=Electron,FinalFlav=Electron] ≈ 0.834 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Electron,FinalFlav=Muon]     ≈ 0.110 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Electron,FinalFlav=Tau]      ≈ 0.055 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Muon,FinalFlav=Electron]     ≈ 0.079 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Muon,FinalFlav=Muon]         ≈ 0.054 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Muon,FinalFlav=Tau]          ≈ 0.867 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Tau,FinalFlav=Electron]      ≈ 0.083 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Tau,FinalFlav=Muon]          ≈ 0.836 atol=0.01
@test test_values[Energy=1,Path=1,InitFlav=Tau,FinalFlav=Tau]           ≈ 0.077 atol=0.01

test_values = Pνν(osc, E, p, anti=true)
@test test_values[1,1,1,1] ≈ 0.976 atol=0.01
@test test_values[1,1,1,2] ≈ 0.013 atol=0.01
@test test_values[1,1,1,3] ≈ 0.011 atol=0.01
@test test_values[1,1,2,1] ≈ 0.012 atol=0.01
@test test_values[1,1,2,2] ≈ 0.663 atol=0.01
@test test_values[1,1,2,3] ≈ 0.325 atol=0.01
@test test_values[1,1,3,1] ≈ 0.013 atol=0.01
@test test_values[1,1,3,2] ≈ 0.324 atol=0.01
@test test_values[1,1,3,3] ≈ 0.664 atol=0.01

h5open("data/refdata.h5", "r") do file
    # Nu-Fit v5.0 Values
    osc_nh = OscillationParameters(3);
    mixingangle!(osc_nh, 1=>2, 5.836e-1);
    mixingangle!(osc_nh, 1=>3, 1.496e-1);
    mixingangle!(osc_nh, 2=>3, 8.587e-1);
    cpphase!(osc_nh, 1=>3, 197 * π / 180);
    masssquareddiff!(osc_nh, 2=>3, -2.517e-3);
    masssquareddiff!(osc_nh, 1=>2, -7.42e-5);

    density = collect(0.1:0.1:40)
    baseline = 10000
    energy = 10 .^ (range(0;stop=2,length=500))
    paths = [Neurthino.Path(d, baseline) for d in density]

    U_nh = PMNSMatrix(osc_nh)
    H_nh = Hamiltonian(osc_nh)

    data_matter_nh = Neurthino.Pνν(U_nh, H_nh, energy, paths; anti=false);

    refdata = read(file, "matter/prob_nh")
    @test data_matter_nh ≈ refdata atol=0.01
end

osc = OscillationParameters(4);

setθ!(osc, 1=>2, 0.59);
setθ!(osc, 1=>3, 0.15);
setθ!(osc, 2=>3, 0.84);
setδ!(osc, 1=>3, 3.86);

setΔm²!(osc, 2=>3, -2.523e-3);
setΔm²!(osc, 1=>2, -7.39e-5);

setθ!(osc, 1=>4, 0);
setθ!(osc, 2=>4, 0.1);
setθ!(osc, 3=>4, 0.);
setΔm²!(osc, 1=>4, -1);

p = Neurthino.Path([13],[12742])
E = [1000, 2000, 3000]
test_values = Pνν(osc, E, p, anti=false);
test_values_anti = Pνν(osc, E, p, anti=true);

@test test_values[1, 1, 2, 2] ≈ 0.984 atol=0.01
@test test_values[2, 1, 2, 2] ≈ 0.999 atol=0.01
@test test_values[3, 1, 2, 2] ≈ 0.997 atol=0.01

@test test_values_anti[1, 1, 2, 2] ≈ 0.896 atol=0.01
@test test_values_anti[2, 1, 2, 2] ≈ 0.006 atol=0.01
@test test_values_anti[3, 1, 2, 2] ≈ 0.995 atol=0.01
