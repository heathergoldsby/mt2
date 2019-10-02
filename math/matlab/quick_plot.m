function d = quick_plot(xstr, ystr, ql, varargin)
%QUICK_PLOT    Intelligently load and plot the specified data files.
%
%    D = QUICK_PLOT(XSTR,YSTR,QL) plots the XSTR and YSTR values from the data 
%    the data files given in QL.  See below for more information on these
%    parameters.
%
%    D = QUICK_PLOT(XSTR,YSTR,QL,YTRANSFORM) plots the XSTR and YSTR values from 
%    the data files given in QL, after first transforming the values with
%    YTRANSFORM.  See below for more information on these parameters.
%
%    D = QUICK_PLOT(XSTR,YSTR,QL,PROP...) plots the XSTR and YSTR
%    values from the data files given in QL.  See below for information on
%    the PROP parameter.
%
%    D: A cell array of the raw (untransformed) data that was used to
%    generate the plot.
%
%    QL: An (m x 2) cell matrix of quick_load-compatible load strings and
%    title strings.  For example:
%      {'"001-structure/ini_.*gen-Genchamp-AvgFit"', 'hn-ini';...
%       '"001-structure/geo1_.*gen-Genchamp-AvgFit"', 'hn-geo1'}
%
%    XSTR, YSTR: Strings that are eval'd on the data returned from
%    quick_load to yield the x and y values of the data to be plotted.
%    For example, if XSTR were 'data(1,:,1)', then d.data(1,:,1) would be
%    used as the x-value of data to be plotted.
%
%    YTRANSFORM: A tranformation function that is applied to the data to be
%    plotted.  This can be used to unwind fitness functions, scale data,
%    etc.  It is a two-parameter function of the form f(i,x), where i is
%    the index number in the QL string, and x is the value to be
%    transformed.  For example, @(i,x) (x./1600) would divide all values by
%    1600.
%
%    PROP: Additional properties that can be set via named parameters.
%    Valid properties include:
%      'plotse': Whether to plot standard error; 0=no, 1=yes (default).
%      'combine': Function used to combine data; defaults to @mean.
%      'xtransform': Function used to transform the x-axis; see YTRANSFORM.
%      'ytransform': See YTRANSFORM, above.
%
%    EXAMPLE:
%      d = quick_plot('data(1,:,1)', 'data(:,:,2)',...
%          {'"001-structure/ini_.*gen-Genchamp-AvgFit"', 'hn-ini';...
%           '"001-structure/geo1_.*gen-Genchamp-AvgFit"', 'hn-geo1';...
%           '"001-structure/3d_.*gen-Genchamp-AvgFit"', '3d';...
%           '"001-structure/mhn_.*gen-Genchamp-AvgFit"', 'ma-hn';...
%           '"001-structure/n_.*gen-Genchamp-AvgFit"', 'n'},...
%          @(i,x) (x./1600));
%      axis([0 200 0 1.05]);
%      xlabel('Generations')
%      ylabel('Fraction of max. fitness');
%      title('Varied neuroevolutionary approaches');
%      quick_export('../doc/gecco2010/figures/001-structure.eps');

% input validation
%
if(nargin < 3)
    error('Wrong number of input arguments: xstr, ystr, and ql are all required');
end

% default values
%
bounds = 'ci'; % show standard error
xcombinefcn = @(x)(x(1,:));
ycombinefcn = @mean; % calculate the mean
xtransform = @(i,x)x; % identity
ytransform = @(i,y)y; % identity
style = {'b','^',':';...
    'r','o',':';...
    'g','*',':';...
    'm','+',':';...
    'k','s',':';...
    'y','v',':'}; % style for lines; color, marker, std err style
markers = true; % show markers?
markers_per_line=12; % markers per line
marker_size=8; % marker size
line_width=1.0; % line width
d = 0.0; % data to be plotted
d_loaded = false; % do we need to load the data?
bootfold = 200; % number of times to sample for bootstrapping
alpha = 0.05; % ci level
show_legend = true; % whether to show the legend
independent = false; % whether to plot series independently
trimend = false; % whether to trim trailing zeros from the data before plotting
% handle varargs
%
if(nargin > 4); 
    iend = size(varargin,2);
    i = 1;
    while(i<=iend)
        if(strcmp(varargin{i}, 'bounds'))
            bounds = varargin{i+1};
            i=i+1;
            if(strcmp(bounds, 'se'))
            elseif(strcmp(bounds, 'ci'))
            elseif(strcmp(bounds, 'none'))
            else
                error(['Unknown bounds type: ' bounds]);
            end
        elseif(strcmp(varargin{i}, 'xtransform'))
            xtransform = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'ytransform'))
            ytransform = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'xcombinefcn'))
            xcombinefcn = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'ycombinefcn'))
            ycombinefcn = varargin{i+1};
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
        elseif(strcmp(varargin{i}, 'bootfold'))
            bootfold = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'alpha'))
            alpha = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'show_legend'))
            show_legend = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'markers'))
            markers = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'independent'))
            independent = varargin{i+1};
            i=i+1;
        elseif(strcmp(varargin{i}, 'trimend'))
            trimend = varargin{i+1};
            i=i+1;
        else
            error(['Unknown parameter "' varargin{i} '" at position ' num2str(i)]);
        end
        i=i+1;
    end
end

% the 'independent' option implies the following options:
if independent
    trimend = true;
    markers = false;
    bounds = 'none';
    show_legend = false;
end

% if the data wasn't passed in, load it:
if(~d_loaded || isempty(d))
    d = quick_multiload(ql);
end

% plot it!
lh = zeros(size(d,1),1);
for i=1:size(d,1)
    x = eval(['d{' num2str(i) '}.' xstr]);
    y = eval(['d{' num2str(i) '}.' ystr]);
    skip = floor(size(x,2)/markers_per_line);

    if independent   
        if size(x,1) ~= size(y,1)
            error('(quick_plot) sizes of x and y in independent data must be equal');
        end
        for j=1:size(y,1)
            yj = y(j,:);
            if trimend
                yj = trim_end(yj);
            end
            xj = x(j,1:length(yj));
            plot(xj, yj, 'linewidth', line_width);
        end
    else
        x_t = xtransform(i,x);
        x = xcombinefcn(x_t);

        y_t = ytransform(i,y);
        % if the number of rows in y is greater than 1, don't transform,
        % as the default behavior for functions like mean will produce a
        % single value, and thus confuse everything.
        if size(y,1) > 1
            y = ycombinefcn(y_t);
        else
            y = y_t;
        end
        p = plot(x, y, style{i,1}, 'linewidth', line_width);
        lh(i) = p(1);
    end
        
    if markers
        plot(x(1:skip:end), y(1:skip:end), [style{i,1} style{i,2}],...
        'markersize', marker_size);
    end
    offset = 1;%floor(skip/2);
    
    % don't try to build error bars if we're only plotting a single row, or
    % if we're plotting independent rows
    if ~strcmp(bounds,'none') && size(y_t,1) > 1
        if(strcmp(bounds, 'se'))
            yse = stderr(y_t(:,offset:skip:end));
            errorbar(x(offset:skip:end), y(offset:skip:end), yse, yse, [style{i,1} style{i,2}]);
        elseif(strcmp(bounds, 'ci'))
            [lb,ub] = all_bootci(y_t(:,offset:skip:end), bootfold, ycombinefcn, alpha);
            ub = ub - y(offset:skip:end)';
            lb = lb - y(offset:skip:end)';
            errorbar(x(offset:skip:end), y(offset:skip:end), lb, ub, [style{i,1} '.'], 'markersize', 1);
        end
    end
end

if(show_legend)
%    [lh, oh, ph, s] = legend(lh, ql{:,2});%,
%    'orientation','horizontal','location','southoutside');
    [lh, oh, ph, s] = legend(lh, ql{:,2},'location','east');
    for i=1:size(d,1)
        set(oh(size(d,1)+i*2),'marker',style{i,2});
    end
end