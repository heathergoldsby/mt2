function H = joint_entropy(X)
%JOINT_ENTROPY    Calculate joint information entropy in bits.
%
%   H = JOINT_ENTROPY(X) calculates the joint entropy H(n0,n1,...) in
%   bits.  The columns of matrix X are the random variables, while the rows
%   of X are the events.  While H(X,Y) is defined as:
%      H(X,Y) = -sum_x sum_y P(x,y) log P(x,y)
%   because we define each row in X as a single event, we are able to 
%   iterate over rows instead of doing the sum_x sum_y.
if ~ismatrix(X)
    error('X must be a matrix.')
end

P = probability_mass(X);
H = -sum(prod([P; log2(P)]));
