import Neurthino: PREM

@test PREM(0) ≈ 13.0885
@test PREM(1000) ≈ 12.87075725
@test PREM(1234) ≈ 12.15988934
@test PREM(6370) ≈ 1.02
