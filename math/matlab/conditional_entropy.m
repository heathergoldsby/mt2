function H = conditional_entropy(X,Y)
%CONDITIONAL_ENTROPY    Calculate conditional information entropy.
%
%   H = CONDITIONAL_ENTROPY(X,Y) calculates the conditional entropy H(X|Y)
%   in bits.  This simply uses the relationship:
%       H(X|Y) = H(X,Y) - H(Y)
H = joint_entropy([X Y]) - entropy(Y);