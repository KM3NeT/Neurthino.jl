using Neurthino

@test Neurthino.cpphases(0) == 0
@test Neurthino.cpphases(1) == 0
@test Neurthino.cpphases(2) == 0
@test Neurthino.cpphases(3) == 1
@test Neurthino.cpphases(4) == 3

osc = Neurthino.OscillationParameters(3)
mixingangle!(osc, 1=>2, 0.59)
Neurthino.mixingangle!(osc, 1=>3, 0.15)
Neurthino.mixingangle!(osc, 2=>3, 0.84)
Neurthino.masssquareddiff!(osc, 2=>3, -2.523e-3)
Neurthino.masssquareddiff!(osc, 1=>2,-7.39e-5)
U = Neurthino.PMNSMatrix(osc)
H = Neurthino.Hamiltonian(osc)

test_values = Neurthino.oscprob(U, H, 1, 1.6e4)
@test test_values[1,1] ≈ 0.142 atol=0.01 
@test test_values[1,2] ≈ 0.420 atol=0.01 
@test test_values[1,3] ≈ 0.438 atol=0.01 
@test test_values[2,2] ≈ 0.256 atol=0.01 
@test test_values[2,3] ≈ 0.323 atol=0.01 
@test test_values[3,3] ≈ 0.239 atol=0.01 

Neurthino.masssquareddiff!(osc, 3=>2, 2.523e-3)
Neurthino.masssquareddiff!(osc, 2=>1, 7.39e-5)

test_values = Neurthino.oscprob(osc, 1, 1.6e4)
@test test_values[1,1] ≈ 0.142 atol=0.01 
@test test_values[1,2] ≈ 0.420 atol=0.01 
@test test_values[1,3] ≈ 0.438 atol=0.01 
@test test_values[2,2] ≈ 0.256 atol=0.01 
@test test_values[2,3] ≈ 0.323 atol=0.01 
@test test_values[3,3] ≈ 0.239 atol=0.01 

Neurthino.cpphase!(osc, 1=>3, 3.86)

U = Neurthino.PMNSMatrix(osc)
H = Neurthino.Hamiltonian(osc)

U_mat, H_mat = Neurthino.MatterOscillationMatrices(U, H, 13.0, 1.0)

matter_test_values = Neurthino.oscprob(U_mat, H_mat, 1, 1.0e4)
@test test_values[1,1] ≈ 0.142 atol=0.01 
@test test_values[1,2] ≈ 0.420 atol=0.01 
@test test_values[1,3] ≈ 0.438 atol=0.01 
@test test_values[2,2] ≈ 0.256 atol=0.01 
@test test_values[2,3] ≈ 0.323 atol=0.01 
@test test_values[3,3] ≈ 0.239 atol=0.01 

U_mat, H_mat = Neurthino.MatterOscillationMatrices(osc, 13.0, 1.0)

matter_test_values = Neurthino.oscprob(U_mat, H_mat, 1, 1.0e4)
@test test_values[1,1] ≈ 0.142 atol=0.01 
@test test_values[1,2] ≈ 0.420 atol=0.01 
@test test_values[1,3] ≈ 0.438 atol=0.01 
@test test_values[2,2] ≈ 0.256 atol=0.01 
@test test_values[2,3] ≈ 0.323 atol=0.01 
@test test_values[3,3] ≈ 0.239 atol=0.01 

@test_logs (:warn, "Mass squared difference fields (partially) overdetermined!") Neurthino.masssquareddiff!(osc, 3=>1, 1)

for i in 3:100
    osc_params_dims = Neurthino.OscillationParameters(i)
    @test_throws ErrorException("Mass squared difference with equal index cannot be modified.") setΔm²!(osc_params_dims, 1=>1, 1) 
    for j in 1:i-1
        @test_throws ErrorException("Mass squared differences not fully determined!") Hamiltonian(osc_params_dims)
        Neurthino.masssquareddiff!(osc_params_dims, j=>j+1, 1.0)
    @test_logs (:warn, "Mass squared difference fields (partially) overdetermined!") Neurthino.masssquareddiff!(osc, 3=>1, 1)
    end
end

