function I = conditional_mutual_information(X,Y,Z)
%CONDITIONAL_MUTUAL_INFORMATION    Calculate conditional mutual information
%   in bits.
%
%   I = CONDITIONAL_MUTUAL_INFORMATION(X,Y,Z) calculates the information 
%   I(X;Y|Z) in bits.  This is defined as:
%       I(X;Y|Z) = H(X,Z) + H(Y,Z) - H(X,Y,Z) - H(Z)
I = joint_entropy([X Z]) + joint_entropy([Y Z]) - joint_entropy([X Y Z]) - entropy(Z);