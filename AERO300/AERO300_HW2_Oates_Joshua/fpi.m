function [x, count] = fpi(g, x_0, TOL)
y_0 = g(x_0);
count = 0;
max_it = 10000;
while (abs(x_0 - y_0) > TOL)
    x_0 = y_0;
    y_0 = g(x_0);
    count = count + 1;
    
    if (count > max_it)
        display('Max number of iterations reached.  Algorithm does not appear to beconverging')
        break;
    end
end
x = x_0;