function [H D C] = ap_joint_entropy(X)
%AP_JOINT_ENTROPY    Calculates joint information entropy among all
%   columns of X.
%
%   [H D C] = ALL_PAIRS_JOINT_ENTROPY(X) H, D, and C are all matrices of
%   size n x n, where n is the number of columns in X.  H(i,j) is the joint
%   entropy between X(:,i) and X(:,j), D is the difference between the
%   joint entropy and the sum of the independent entropies, and C is true
%   if i and j are conditional on each other.  H, D, and C are upper
%   triangular and the diagonal is not computed.

cols = size(X,2);
H = zeros(cols);
D = zeros(cols);
C = zeros(cols);

for i=1:cols
    for j=(i+1):cols
        [H(i,j) D(i,j) C(i,j)] = joint_entropy(X(:,i), X(:,j));
    end
end
