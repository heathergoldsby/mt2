function v = deg2vec(theta)
%DEG2VEC    Translates degrees to a unit vector.
%
%    V = DEG2VEC(THETA) returns a unit vector that is theta degrees from
%    the vector (1,0).
v = [cosd(theta) sind(theta)];