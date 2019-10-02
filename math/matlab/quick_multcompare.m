function quick_multcompare(ystr, ql, varargin)
%QUICK_MULTCOMPARE    Intelligently load and perform a multiple comparison
%    the specified data files.
%

% input validation
%
if(nargin < 2)
    error('Wrong number of input arguments: xstr, ystr, and ql are all required');
end

% default values
%
d = 0.0; % data to be plotted
d_loaded = false; % do we need to load the data?
alpha = 0.05; % alpha to test against

% handle varargs
%
if(nargin > 4); 
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'data'))
            d = varargin{i+1};
            d_loaded = true;
            i=i+1;
        elseif(strcmp(varargin{i}, 'alpha'))
            alpha = varargin{i+1};
            i=i+1;
        else
            error(strcat('Unknown parameter at position ', num2str(i)));
        end
        i=i+1;
    end
end    

% if the data wasn't passed in, load it:
if(~d_loaded || isempty(d))
    d = quick_multiload(ql);
end

% do the multiple comparison!
%
n = size(d,1); % number of treatments
m = size(eval(['d{1}.' ystr]),1); % number of data points within treatments

mc = zeros(m,n);

% each column is data from one treatment; 
for i=1:n
    mc(:,i) = eval(['d{' num2str(i) '}.' ystr]);
end

[p, table, stats] = kruskalwallis(mc);
multcompare(stats,'alpha', alpha);