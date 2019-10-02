function [P V C] = probability_mass(X)
%PROBABILIITY_MASS    Calculates the probability mass of events in X.
%
%   [P V C] = PROBABILIITY_MASS(X) calculates the set of unique values (V), 
%   their counts (C), and their probability (P) in the random variable X.
%   Each row of X is treated as an event.  In all cases, V contains a 
%   character representation of the event, rather than the event as-is from
%   X.
pmf = containers.Map('keytype','char','valuetype','double');
events = size(X,1);

for i = 1:events
    k = num2str(X(i,:));
    if isKey(pmf,k)
        pmf(k) = pmf(k) + 1;
    else
        pmf(k) = 1;
    end
end

V = keys(pmf);
C = zeros(1,length(V));
P = zeros(1,length(V));
for i=1:length(V)
    C(i) = pmf(V{i});
    P(i) = C(i)/events;
end