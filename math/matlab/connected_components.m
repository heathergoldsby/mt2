function [cc, nc] = connected_components(A)
% A == adjacency matrix, cc == component labeling

spA = spones(A) + speye(size(A));
[p,q,r,s] = dmperm(spA);

cc = zeros(1,size(A,1));
for j=1:(length(r)-1)
    cc(p(r(j):(r(j+1)-1))) = r(j);
end

nc = length(r);