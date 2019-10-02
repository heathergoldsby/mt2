function d = quick_ciplot(xstr, ystr, ql, varargin)
%QUICK_CIPLOT    Intelligently load and plot the specified data files, with
%    confidence intervals.
%

% input validation
%
if(nargin < 3)
    error('Wrong number of input arguments: xstr, ystr, and ql are all required');
end

% default values
%
xtransform = @(x)x; % identity
ytransform = @(i,x)x; % identity
style = {'b','^',':';...
    'r','o',':';...
    'g','*',':';...
    'm','+',':';...
    'k','s',':';...
    'y','v',':'}; % style for lines; color, marker, std err style
d = 0.0; % data to be plotted
d_loaded = false; % do we need to load the data?
skip = 4; % number of points to skip between bootstraps
bootfold = 200; % number of bootstraps to run
bootfcn = @mean; % function to calc bootstrap on

% handle varargs
%
if(nargin > 4); 
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'xtransform'))
            xtransform = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'ytransform'))
            ytransform = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'style'))
            style = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'data'))
            d = varargin{i+1};
            d_loaded = true;
            i=i+1;
        elseif(strcmp(varargin{i}, 'skip'))
            skip = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'bootfold'))
            bootfold = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'bootfcn'))
            bootfcn = varargin{i+1};
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

% plot it!
newfigure();
lh = zeros(size(d,1),1);
for i=1:size(d,1)
    x = eval(['d{' num2str(i) '}.' xstr]);
    x = xtransform(x);    
    
    y = eval(['d{' num2str(i) '}.' ystr]);
    y = ytransform(i,y);
    
    p = plot(x,bootfcn(y), style{i,1}, 'linewidth', 1.5);
    lh(i) = p(1);
    
    [lb,ub] = all_bootci(y(:,1:skip:end), bootfold, bootfcn, 0.05);
    ciplot(lb, ub, x(1:skip:end), style{i,1});
end

legend(lh, ql{:,2}, 'orientation','horizontal','location','southoutside');
