# part 1
class vehicle():
    biome = "Universe"
    def __init__(self,make,model,year):
        self.make = make
        self.model = model
        self.year = year
    def __str__(self):
        return(f"{self.year} \n{self.make} \n{self.model}")


class car(vehicle):
    numWheels = 4
    biome = "Road"
    title = "clean"
    def __init__(self,make,model,year,color):
        super().__init__(make,model,year)
        self.color = color
    def __str__(self):
        return(f"{self.color} \n{self.year} \n{self.make} \n{self.model}")
    def __add__(self,other):
        self.title = "totaled"
        other.title = "totaled"


class plane(vehicle):
    biome = "Sky"
    def __init__(self,make,model,year,numPassengers):
        super().__init__(make,model,year)
        self.numPassengers = numPassengers
    def airMarshal(self):
        self.numPassengers -= 1




