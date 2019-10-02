function y = scale(x, m)
%SCALE    Scales a vector to the given magnitude.
%
%    Y = SCALE(X) returns the vector X scaled to magnitude 1.0.  This is
%    equivalent to norm(X).
%
%    Y = SCALE(X, M) returns the vector X scaled to magnitude M.
if nargin < 2, m=1.0; end
y = x ./ norm(x) .*m;