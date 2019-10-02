function [lb, ub] = all_bootci(y, nboot, bootfun, alpha)
% ALL_CI
%
% y must be a matrix...
assert(ndims(y) == 2);

lb = zeros(size(y,2),1);
ub = zeros(size(y,2),1);

for i = 1:size(y,2)
    if mean(y(:,i)) ~= -Inf
        ci = bootci(nboot, {bootfun, y(:,i)}, 'alpha', alpha);
        lb(i) = ci(1);
        ub(i) = ci(2);
    end
end