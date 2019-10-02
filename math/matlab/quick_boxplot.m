function [d, l] = quick_boxplot(ystr, ql, varargin)
%QUICK_BOXPLOT    Intelligently load and boxplot the specified data files.

% input validation
if(nargin < 2)
    error('Wrong number of input arguments');
end

% defaults
ytransform = @(i,x)x;
d = cell(size(ql,1),1);
d_loaded = false;

% handle varargs:
if(nargin > 2); 
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'ytransform'))
            ytransform = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'data'))
            d = varargin{i+1};
            d_loaded = true;
            i=i+1;
        else
            error(['Unknown parameter at position ' num2str(i)]);
        end
        i=i+1;
    end
end

% have we loaded data?
if(~d_loaded)
    for i=1:size(ql,1)
        d{i} = quick_load(ql{i,1});
    end
end

% after all that, build the boxplot:
exa = eval(['d{1}.' ystr]);
m = zeros(size(exa,1),size(d,1));

l = cell(size(ql,1),1);
for i=1:size(d,1)
    y = eval(['d{' num2str(i) '}.' ystr]);
    y = ytransform(i,y);
    m(1:size(y,1),i) = y;
    l(i) = ql(i,2);
end
%display(l);
[fig axes] = newfigure();
boxplot(m, 'notch', 'on', 'labels', l);
box('off');