# part 1
import Joshua_Oates_ICA4 as j

myv = j.vehicle("NASA","Saturn-5",1969)
print(myv)
print(myv.biome)
print("\n")


myc = j.car("Honda","Accord",2004,"green")
print(myc)
print(myc.biome)
print("")
yourc = j.car("Mazda","Protege",1995,"white")
print(yourc)
print(yourc.title)
myc+yourc
print(yourc.title)
print("\n")



myp = j.plane("Cessna","152",1995,4)
print(myp)
print(myp.numPassengers)
print(myp.biome)
print("")
myp.airMarshal()
print(myp.numPassengers)
print("\n")




