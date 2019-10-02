function qh = hintze_modularity(A, Ts, c)
% A == adjacency matrix, Ts == partitions, c == number of components

qh = zeros(size(Ts,1),1);
vertices = size(A,1); % num vertices

for p=1:size(Ts,1) % partitions
    m = sum(sum(A))/2; % sum of all edges
    nc = c(p); % number of components
    
    s = 0.0;    
    for i=1:vertices
        for j=(i+1):vertices
            if Ts(p,i) == Ts(p,j)
                s = s + A(i,j);
            else
                s = s - A(i,j) / (nc-1);
            end
        end
    end
    qh(p) = s / m;
end
