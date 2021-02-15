using Neurthino

osc = OscillationParameters(3)
setθ!(osc, 1=>2, 0.59);
setθ!(osc, 1=>3, 0.15);
setθ!(osc, 2=>3, 0.84);
setδ!(osc, 1=>3, 3.86);
setΔm²!(osc, 2=>3, -2.523e-3);
setΔm²!(osc, 1=>2, -7.39e-5);

U_mat, H_mat = MatterOscillationMatrices(osc, 1, 13);
test_values = Pνν(U_mat, H_mat, 1, 10000)
@test test_values[1,1] ≈ 0.834 atol=0.01                                        
@test test_values[1,2] ≈ 0.083 atol=0.01                                        
@test test_values[1,3] ≈ 0.083 atol=0.01                                        
@test test_values[2,1] ≈ 0.108 atol=0.01                                        
@test test_values[2,2] ≈ 0.052 atol=0.01                                        
@test test_values[2,3] ≈ 0.840 atol=0.01                                        
@test test_values[3,1] ≈ 0.058 atol=0.01                                        
@test test_values[3,2] ≈ 0.865 atol=0.01                                        
@test test_values[3,3] ≈ 0.077 atol=0.01  

p = Neurthino.Path([13],[10000])
E = [1,2]
test_values = Pνν(osc, E, p)
@test test_values[1,1,1,1] ≈ 0.834 atol=0.01                                        
@test test_values[1,1,1,2] ≈ 0.083 atol=0.01                                        
@test test_values[1,1,1,3] ≈ 0.083 atol=0.01                                        
@test test_values[1,1,2,1] ≈ 0.108 atol=0.01                                        
@test test_values[1,1,2,2] ≈ 0.052 atol=0.01                                        
@test test_values[1,1,2,3] ≈ 0.840 atol=0.01                                        
@test test_values[1,1,3,1] ≈ 0.058 atol=0.01                                        
@test test_values[1,1,3,2] ≈ 0.865 atol=0.01                                        
@test test_values[1,1,3,3] ≈ 0.077 atol=0.01  


p = Neurthino.Path([13,13,13,13],[2500,2500,2500,2500])
E = 1

test_values = Pνν(osc, E, p)
@test test_values[1,1,1,1] ≈ 0.834 atol=0.01                                        
@test test_values[1,1,1,2] ≈ 0.083 atol=0.01                                        
@test test_values[1,1,1,3] ≈ 0.083 atol=0.01                                        
@test test_values[1,1,2,1] ≈ 0.108 atol=0.01                                        
@test test_values[1,1,2,2] ≈ 0.052 atol=0.01                                        
@test test_values[1,1,2,3] ≈ 0.840 atol=0.01                                        
@test test_values[1,1,3,1] ≈ 0.058 atol=0.01                                        
@test test_values[1,1,3,2] ≈ 0.865 atol=0.01                                        
@test test_values[1,1,3,3] ≈ 0.077 atol=0.01  


h5open("data/refdata.h5", "r") do file
    # Nu-Fit v5.0 Values
    osc_nh = OscillationParameters(3);
    mixingangle!(osc_nh, 1=>2, 5.836e-1);
    mixingangle!(osc_nh, 1=>3, 1.496e-1);
    mixingangle!(osc_nh, 2=>3, 8.587e-1);
    cpphase!(osc_nh, 1=>3, 197 * π / 180);
    masssquareddiff!(osc_nh, 2=>3, -2.517e-3);
    masssquareddiff!(osc_nh, 1=>2, -7.42e-5);

    density = collect(0.1:1.0:40)
    baseline = 10000
    energy = 10 .^ (range(0;stop=2,length=50))
    paths = [Neurthino.Path(d, baseline) for d in density]
    
    U_nh = PMNSMatrix(osc_nh)
    H_nh = Hamiltonian(osc_nh)

    data_matter_nh = Neurthino.oscprob(U_nh, H_nh, energy, paths; anti=false);

    refdata = read(file, "matter/prob_nh")
    @test data_matter_nh ≈ refdata atol=0.01
end
