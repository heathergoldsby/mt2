function I = ap_mutual_information(X)
%AP_MUTUAL_INFORMATION    Calculates mutual information among all pairs of 
%   columns of X.
%
%   I = AP_MUTUAL_INFORMATION(X) calculates mutual information for all
%   pairs of columns of X.  I is an n x n matrix, where n is the number of
%   columns in X.  I is upper triangular, and the diagonal is not computed.
cols = size(X,2);
I = zeros(cols);

for i=1:cols
    for j=(i+1):cols
        I(i,j) = mutual_information(X(:,i), X(:,j));
    end
end
