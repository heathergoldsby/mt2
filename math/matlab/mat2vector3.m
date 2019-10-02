function [x y z] = mat2vector3(X)

r = rows(X); c = cols(X);
vsize = r*c;
x = zeros(vsize,1);
y = zeros(vsize,1);
z = zeros(vsize,1);

for i=1:r
    for j=1:c
        idx = (i-1)*c + j;
        x(idx) = j;
        y(idx) = i;
        z(idx) = X(i,j);
    end
end