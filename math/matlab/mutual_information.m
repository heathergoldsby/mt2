function I = mutual_information(X,Y)
%MUTUAL_INFORMATION    Calculate mutual information in bits.
%
%   I = MUTUAL_INFORMATION(X,Y) calculates the information I(X;Y) in bits.
%   This is defined as:
%       I(X;Y) = H(X) + H(Y) - H(X,Y)
I = entropy(X) + entropy(Y) - joint_entropy([X Y]);