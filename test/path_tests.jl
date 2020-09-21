import Neurthino

@test Neurthino.tracklength(60u"°", Neurthino.EARTH_RADIUS / 2) ≈ 4149.99u"km" atol=0.01u"km"
@test Neurthino.tracklength(40u"°", Neurthino.EARTH_RADIUS / 2) ≈ 3592.76u"km" atol=0.01u"km"
@test Neurthino.tracklength(15u"°", Neurthino.EARTH_RADIUS / 2) ≈ 3240.47u"km" atol=0.01u"km"
@test Neurthino.tracklength(110u"°", Neurthino.EARTH_RADIUS / 2) ≈ 6713.49u"km" atol=0.01u"km" 
@test Neurthino.tracklength(140u"°", Neurthino.EARTH_RADIUS / 2) ≈ 8473.23u"km" atol=0.01u"km" 

@test Neurthino.tracklength(150u"°", 1860u"km") ≈ 9865.03u"km" atol=0.01u"km" 
@test Neurthino.tracklength(170u"°", 1860u"km") ≈ 10765.13u"km" atol=0.01u"km" 
@test Neurthino.tracklength(130u"°", 5600u"km") ≈ 6839.15u"km" atol=0.01u"km" 



