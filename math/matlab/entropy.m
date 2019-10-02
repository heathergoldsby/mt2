function H = entropy(X)
%ENTROPY    Calculate information entropy in bits.
%
%   H = ENTROPY(X) calculates the entropy H(X) in bits.  If X is a row or column 
%   vector, H is a scalar value.  If X is a matrix, then entropy is
%   calculated for each column and H is a row vector. In all cases, X is 
%   assumed to be a vector/matrix of events -- That is, the value of any 
%   element in X is immaterial, and only the relative frequencies of values
%   are used to calculate entropy.
if isvector(X)
    H = vector_entropy(X);
elseif ~ismatrix(X)
    error('X must be either a vector or matrix of events.');
else
    H = zeros(1,size(X,2));
    for i=1:length(H)
        H(i) = vector_entropy(X(:,i));
    end
end
end
        
function h = vector_entropy(x)
    P = probability_mass(x);
    h = -sum(prod([P; log2(P)]));
end