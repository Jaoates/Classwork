A = [1];
B = [1,2];
C = [[1,2,3];[4,5,6]];

n=1;
for i= 1:4
    for j = 1:4
        for k = 1:4
            D(i,j,k)=n;
            n=n+1;
        end
    end
end

n=1;
for i= 1:5
    for j = 1:5
        for k = 1:5
           for t=1:5
            E(i,j,k,t)=n;
            n=n+1;
           end
        end
    end
end