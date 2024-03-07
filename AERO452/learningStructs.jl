struct Cat
    name::String
    color::String
    age::Number
    sound::String
end

function Cat( # there are too many args to have on one line
    name::String;
    color = "Plain",
    age = NaN,
    sound = "Meow"
    )
    Cat(name,color,age,sound)
end

function Cat(s1::Cat,s2::Cat)
    name = s1.name*"-"*s2.name
    color = s1.color*" and "*s2.color
    age = 0
    Cat(name,color = color,age = age)
end

function Base.:+(x::Cat,y::Cat)
    Cat(x,y)
end

function Base.show(io::IO,s::Cat) # ive overwritten the show function for Cat types
     print("The ",s.color," cat named ",s.name)
end

function makeSound(s::Cat)
    println(s.sound)
end

s = Cat("Soph","Orange",22,"Sup?")
j = Cat("Josh",sound = "Hi there",)

println(s)
println(j)

makeSound(s)
makeSound(j)
sj = s+j
println(sj)