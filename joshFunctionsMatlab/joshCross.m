function mx = joshCross(m)
% takes a column vector and returns the associated 'cross' matrix such that
% mx*b == cross(m,b)
arguments
    m (3,1) 
end
if isa(m,"double") % overloaded to handle either symbolic or double type vectors
    mx = zeros(3);
elseif isa(m,"sym")
    syms mx [3 3]
else
    throw(MException("joshCross:invalidInput","m must be type sym or double"))
end
    for i = 1:3
        mx(i,i) = 0;
    end
    mx(1,2) = -m(3);
    mx(1,3) =  m(2);
    mx(2,3) = -m(1);

    mx(2,1) =  m(3);
    mx(3,1) = -m(2);
    mx(3,2) =  m(1);
end