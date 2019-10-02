function Qn = newman_modularity(A, Ts)
% A == adjacency matrix, Ts == partitions

Qn = zeros(size(Ts,1),1);
vertices = size(A,1); % num vertices

for p=1:size(Ts,1) % partitions
    m = sum(sum(A))/2; % sum of all edges
    s = 0.0;
    for i=1:vertices
        for j=(i+1):vertices
            if Ts(p,i) == Ts(p,j)
                ki = sum(A(:,i))/2;
                kj = sum(A(:,j))/2;
                aij = A(i,j);
                s = s + aij - (ki*kj)/(2*m);
            end
        end
    end
    Qn(p) = 1.0 / (4.0 * m) * s;
end
