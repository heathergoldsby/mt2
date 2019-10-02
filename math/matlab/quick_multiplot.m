function d = quick_multiplot(xstr, ystr, ql, varargin)
%QUICK_MULTIPLOT    Intelligently load and plot the specified data files.
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
combinefcn = @mean; % calculate the mean
style = {'b','^',':';...
    'r','o',':';...
    'g','*',':';...
    'm','+',':';...
    'k','s',':';...
    'y','v',':'}; % style for lines; color, marker, std err style
markers_per_line=12; % markers per line
marker_size=8; % marker size
line_width=1.0; % line width
d = 0.0; % data to be plotted
d_loaded = false; % do we need to load the data?

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
        elseif(strcmp(varargin{i}, 'combinefcn'))
            combinefcn = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'style'))
            style = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'markers_per_line'))
            markers_per_line = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'marker_size'))
             marker_size = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'line_width'))
            line_width = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'data'))
            d = varargin{i+1};
            d_loaded = true;
            i=i+1;
        else            
            error(['Unknown parameter "' varargin{i} '" at position ' num2str(i)]);
        end
        i=i+1;
    end
end    

% if the data wasn't passed in, load it:
if(~d_loaded || isempty(d))
    d = quick_multiload(ql);
end

% plot it!
[fig axes] = newfigure();
%lh = [ ];
lh = zeros(size(d,1),1);
for i=1:size(d,1)
    x = eval(['d{' num2str(i) '}.' xstr]);
    x = xtransform(x);    
    
    y = eval(['d{' num2str(i) '}.' ystr]);
    y_t = ytransform(i,y);
    y = combinefcn(y_t);
    skip = floor(size(x,2)/markers_per_line);
    
    p = plot(x, y, style{i,1},...
        x(1:skip:end), y(1:skip:end), [style{i,1} style{i,2}],...
        'markersize', marker_size, 'linewidth', line_width);
    lh(i) = p(1);

    plot(x, y_t, [style{i,1}], 'linewidth', line_width/2);
end

legend(lh, ql{:,2}, 'orientation','horizontal','location','southoutside');
[lh, oh, ph, s] = legend(axes);

for i=1:size(d,1)
    set(oh(size(d,1)+i*2),'marker',style{i,2});
end