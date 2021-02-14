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
