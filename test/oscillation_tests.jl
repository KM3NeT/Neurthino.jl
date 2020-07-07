import Neurthino

@test Neurthino.number_cp_phases(0) == 0
@test Neurthino.number_cp_phases(1) == 0
@test Neurthino.number_cp_phases(2) == 0
@test Neurthino.number_cp_phases(3) == 1
@test Neurthino.number_cp_phases(4) == 3


