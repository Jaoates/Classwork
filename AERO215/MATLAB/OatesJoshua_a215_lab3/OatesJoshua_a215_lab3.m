%% Joshua Oates - a215 - Fall 2021 - Lab3
clear all
close all
clc

%% Part1
%Using “for” loops, create a 3x4 matrix of random numbers. Display the 
%smallest value of each row, and display which element it is (row, column). 

%create some colunms by iterating and adding a colunm to a matrix
for row = 1:3
    smallest = 1;
    colNumber = 0;
    %create each col by iterating and adding random number to it
    for col = 1:4
        random = rand();
        randMat(row,col)=random;
        %check to see if if the current random number is smaller than the
        %previous smallest value or not
        if (random < smallest);
            %save the new smallest value 
            smallest = random;
            %save the column number
            colNumber = col;
        end
    end
    disp("The smallest value in row " + row + " is " + smallest +" in col " + colNumber);
end
disp (randMat);

%% Part 2
%Say you are comparing cars;  one is a grey 2013 Toyota Prius Plug-in 
% hatchback with 134 horsepower that weighs 3165 lbf. The other is a white 
% 2015 Honda Civic LX sedan with 143 horsepower that weighs 2800 lbf. Create 
% one structure for each car, with six elements: color, year, make, model, 
% body style,  horsepower, and weight. Write code using the structures to 
% solve for and report which has greater  horsepower to weight ratio

car1.color = "grey";
car1.year = 2013;
car1.make = "Toyota";
car1.model = "Prius";
car1.bodyStyle = "hatchback";
car1.horsepower = 134;
car1.weight = 3165;
car1.powerToWeight = (car1.horsepower / car1.weight);

car2.color = "white";
car2.year = 2015;
car2.make = "Honda";
car2.model = "Civic";
car2.bodyStyle = "sedan";
car2.horsepower = 143;
car2.weight = 2800;
car2.powerToWeight = (car2.horsepower / car2.weight);

%make a comparison of the two cars' power to weight ratio and print the one
%that has a better ratio.
if car1.powerToWeight > car2.powerToWeight
    disp(car1.make + " " +  car1.model + " " + car1.year + " has a better power to wieght by " + (car1.powerToWeight-car2.powerToWeight) + "hp/lbf")
elseif car1.powerToWeight == car2.powerToWeight
    disp("the cars have the same power to weight")
else
    disp(car2.make + " " +  car2.model + " " + car2.year + " has a better power to wieght by " + (car2.powerToWeight-car1.powerToWeight) + "hp/lbf")
end
