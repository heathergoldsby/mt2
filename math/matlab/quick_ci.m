function [lb,m,ub] = quick_ci(x, varargin)
%QUICK_CI    Quickly calculate a bootstrapped confidence interval for x.

if(nargin < 1)
    error('Wrong number of input arguments: x is required');
end

% default values
%
combinefcn = @mean; % calculate the mean
bootfold = 200; % number of times to sample for bootstrapping
alpha = 0.05; % ci level

% handle varargs
%
if(nargin > 1);
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'combinefcn'))
            combinefcn = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'bootfold'))
            bootfold = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'alpha'))
            alpha = varargin{i+1};
            i=i+1;
        else
            error(['Unknown parameter "' varargin{i} '" at position ' num2str(i)]);
        end
        i=i+1;
    end
end

[lb,ub] = all_bootci(x, bootfold, combinefcn, alpha);
m = combinefcn(x);
