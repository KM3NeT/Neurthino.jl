import Neurthino

@test Neurthino.tracklength(0, 2.0) ≈ 2.0 atol=0.01 
@test Neurthino.tracklength(pi/2, 1.0) ≈ 112.88 atol=0.01

@test Neurthino.tracklength(pi*1/3, Neurthino.EARTH_RADIUS / 2) ≈ 4149.99 atol=0.01
@test Neurthino.tracklength(0.7, Neurthino.EARTH_RADIUS / 2) ≈ 3595.04 atol=0.01
@test Neurthino.tracklength(pi/6, Neurthino.EARTH_RADIUS / 2) ≈ 3409.97 atol=0.01
@test Neurthino.tracklength(1.91, Neurthino.EARTH_RADIUS / 2) ≈ 6678.27 atol=0.01 
@test Neurthino.tracklength(2.44, Neurthino.EARTH_RADIUS / 2) ≈ 8463.26 atol=0.01 

@test Neurthino.tracklength(2.62, 1860) ≈ 9872.5 atol=1.0 
@test Neurthino.tracklength(2.97, 1860) ≈ 10769.02 atol=1.0 
@test Neurthino.tracklength(2.27, 5600) ≈ 6839.15 atol=1.0 

path = Neurthino.prempath(pi, 1.5)

@test path.baseline[1] ≈ 128.69 atol=0.1 

@test path.density[1] ≈ 3.38 atol=0.1
@test path.density[99] ≈ 3.38 atol=0.1
@test path.density[50] ≈ 13.09 atol=0.1
@test path.density[20] ≈ 5.38 atol=0.1
@test path.density[80] ≈ 5.38 atol=0.1
