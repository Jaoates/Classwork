function [f3] = Funtest()
f1 = @(x) x^2;
f2 = @(x) -f1(x);
    function ret = f3(x)
        ret = x^.5
    end
end