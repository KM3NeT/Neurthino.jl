import Neurthino

@test Neurthino.cpphases(0) == 0
@test Neurthino.cpphases(1) == 0
@test Neurthino.cpphases(2) == 0
@test Neurthino.cpphases(3) == 1
@test Neurthino.cpphases(4) == 3

osc = Neurthino.OscillationParameters(3)
Neurthino.mixingangle!(osc, (1,2), 0.59)
Neurthino.mixingangle!(osc, (1,3), 0.15)
Neurthino.mixingangle!(osc, (2,3), 0.84)
Neurthino.masssquareddiff!(osc, (1,3), -2.523e-3)
Neurthino.masssquareddiff!(osc, (2,3), -2.523e-3)
Neurthino.masssquareddiff!(osc, (1,2),-7.39e-5)
U = Neurthino.PMNSMatrix(osc)
H = Neurthino.Hamiltonian(osc)

test_values = Neurthino.transprob(U, H, 1, 1.6e4)
@test test_values[1,1] ≈ 0.360 atol=0.01 
@test test_values[1,2] ≈ 0.437 atol=0.01 
@test test_values[1,3] ≈ 0.202 atol=0.01 
@test test_values[2,2] ≈ 0.330 atol=0.01 
@test test_values[2,3] ≈ 0.233 atol=0.01 
@test test_values[3,3] ≈ 0.565 atol=0.01 

test_values = Neurthino.transprob(osc, 1, 1.6e4)
@test test_values[1,1] ≈ 0.360 atol=0.01 
@test test_values[1,2] ≈ 0.437 atol=0.01 
@test test_values[1,3] ≈ 0.202 atol=0.01 
@test test_values[2,2] ≈ 0.330 atol=0.01 
@test test_values[2,3] ≈ 0.233 atol=0.01 
@test test_values[3,3] ≈ 0.565 atol=0.01 

Neurthino.cpphase!(osc, (1,3), 3.86)

U = Neurthino.PMNSMatrix(osc)
H = Neurthino.Hamiltonian(osc)

H_mat, U_mat = Neurthino.MatterOscillationMatrices(U, H, 13)

matter_test_values = Neurthino.transprob(U_mat, H_mat, 1, 1.0e4)
@test test_values[1,1] ≈ 0.360 atol=0.01 
@test test_values[1,2] ≈ 0.428 atol=0.01 
@test test_values[1,3] ≈ 0.202 atol=0.01 
@test test_values[2,2] ≈ 0.330 atol=0.01 
@test test_values[2,3] ≈ 0.233 atol=0.01 
@test test_values[3,3] ≈ 0.565 atol=0.01 

H_mat, U_mat = Neurthino.MatterOscillationMatrices(osc, 13)

matter_test_values = Neurthino.transprob(U_mat, H_mat, 1, 1.0e4)
@test test_values[1,1] ≈ 0.360 atol=0.01 
@test test_values[1,2] ≈ 0.428 atol=0.01 
@test test_values[1,3] ≈ 0.202 atol=0.01 
@test test_values[2,2] ≈ 0.330 atol=0.01 
@test test_values[2,3] ≈ 0.233 atol=0.01 
@test test_values[3,3] ≈ 0.565 atol=0.01 
