function I = joint_mutual_information(X,Y)
%JOINT_MUTUAL_INFORMATION    Calculates join mutual information in bits 
%   between the joint distribution described by the columns {x1,x2,...xn} 
%   of X and Y.
%
%   I = JOINT_MUTUAL_INFORMATION(X,Y) calculates joint mutual information 
%   in bits.  This is defined as:
%      I(x1,x2,...xn; Y) = H(x1,x2,...xn) + H(Y) - H(x1,x2,...xn,Y)
I = joint_entropy(X) + entropy(Y) - joint_entropy([X Y]);