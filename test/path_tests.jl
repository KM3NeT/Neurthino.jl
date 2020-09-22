import Neurthino

@test Neurthino.tracklength(0u"°", 2.0u"km") ≈ 2.0u"km" atol=0.01u"km" 
@test Neurthino.tracklength(90u"°", 1.0u"km") ≈ 112.88u"km" atol=0.01u"km"

@test Neurthino.tracklength(60u"°", Neurthino.EARTH_RADIUS / 2) ≈ 4149.99u"km" atol=0.01u"km"
@test Neurthino.tracklength(40u"°", Neurthino.EARTH_RADIUS / 2) ≈ 3592.76u"km" atol=0.01u"km"
@test Neurthino.tracklength(15u"°", Neurthino.EARTH_RADIUS / 2) ≈ 3240.47u"km" atol=0.01u"km"
@test Neurthino.tracklength(110u"°", Neurthino.EARTH_RADIUS / 2) ≈ 6713.49u"km" atol=0.01u"km" 
@test Neurthino.tracklength(140u"°", Neurthino.EARTH_RADIUS / 2) ≈ 8473.23u"km" atol=0.01u"km" 

@test Neurthino.tracklength(150u"°", 1860u"km") ≈ 9865.03u"km" atol=0.01u"km" 
@test Neurthino.tracklength(170u"°", 1860u"km") ≈ 10765.13u"km" atol=0.01u"km" 
@test Neurthino.tracklength(130u"°", 5600u"km") ≈ 6839.15u"km" atol=0.01u"km" 

distances, densities = Neurthino.prempath(180u"°", 1500u"m")

for d in distances
    @test distances[1] ≈ 128.69u"km" atol=0.1u"km" 
end

@test densities[1] ≈ 3.38u"g/cm^3" atol=0.1u"g/cm^3"
@test densities[99] ≈ 3.38u"g/cm^3" atol=0.1u"g/cm^3"
@test densities[50] ≈ 13.09u"g/cm^3" atol=0.1u"g/cm^3"
@test densities[20] ≈ 5.38u"g/cm^3" atol=0.1u"g/cm^3"
@test densities[80] ≈ 5.38u"g/cm^3" atol=0.1u"g/cm^3"
