import Neurthino

@test Neurthino.number_cp_phases(0) == 0
@test Neurthino.number_cp_phases(1) == 0
@test Neurthino.number_cp_phases(2) == 0
@test Neurthino.number_cp_phases(3) == 1
@test Neurthino.number_cp_phases(4) == 3

osc = Neurthino.OscillationParameters(3)
osc.mixing_angles[1,2] = 0.59
osc.mixing_angles[1,3] = 0.15
osc.mixing_angles[2,3] = 0.84
osc.mass_squared_diff[1,3] = -2.523e-3
osc.mass_squared_diff[2,3] = -2.523e-3
osc.mass_squared_diff[1,2] = -7.39e-5
U = Neurthino.PMNSMatrix(osc)
H = Neurthino.Hamiltonian(osc)

test_values = Neurthino.transprob(U, H, 1, 1.6e4)
@test test_values[1,1] ≈ 0.360 atol=0.01 
@test test_values[1,2] ≈ 0.437 atol=0.01 
@test test_values[1,3] ≈ 0.202 atol=0.01 
@test test_values[2,2] ≈ 0.330 atol=0.01 
@test test_values[2,3] ≈ 0.233 atol=0.01 
@test test_values[3,3] ≈ 0.565 atol=0.01 
