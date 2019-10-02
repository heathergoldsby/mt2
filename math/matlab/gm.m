function Y = gm(X)
% GM    Safe geometric mean.

Y = zeros(1,size(X,2));
for j=1:size(X,2)
    n = size(X,1);
    x = X(1,j) ^ (1/n);
    for i=2:size(X,1)
%        x = sqrt(x * X(i,j));
        x = x * (X(i,j) ^ (1/n));
    end
    Y(j) = x;
end