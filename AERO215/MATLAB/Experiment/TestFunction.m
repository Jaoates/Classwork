function [outputArg1,outputArg2,out3] = TestFunction(inputArg1,inputArg2,number)

%   how to cheat overload a function. this function can take 2 or 3
%   parameters
    arguments
       inputArg1
       inputArg2
       number = [] %creates a default value of empty var for number if number is not provided. 
    end


outputArg1 = inputArg1;
outputArg2 = inputArg2;
out3 = number;
end

