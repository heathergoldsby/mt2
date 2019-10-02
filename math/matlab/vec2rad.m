function r = vec2rad(x)
%VEC2RAD    Translates a vector to radians.
%
%    R = VEC2RAD(X) converts X to a unit vector uX, and returns the angle
%    (theta) between (1,0) and uX in radians.
r = acos(x(1)/norm(x));