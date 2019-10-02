function d = quick_plotci(xstr, ystr, ql, xform)
%QUICK_PLOTCI    Intelligently load and plot the specified data files, with
%    confidence intervals.
%
%    D = QUICK_PLOT(XSTR,YSTR,QL) plots the XSTR and YSTR values from the data 
%    the data files given in QL.  See below for more information on these
%    parameters.
%
%    D = QUICK_PLOT(XSTR,YSTR,QL,XFORM) plots the XSTR and YSTR values from 
%    the data files given in QL, after first transforming the values with
%    XFORM.  See below for more information on these parameters.
%
%    D: A cell array of the raw (untransformed) data that was used to
%    generate the plot.
%
%    QL: An (m x 2) cell matrix of quick_load-compatible load strings and
%    title strings.  For example:
%      {'001-structure/ini_.*gen-Genchamp-AvgFit', 'hn-ini';...
%       '001-structure/geo1_.*gen-Genchamp-AvgFit', 'hn-geo1'}
%
%    XSTR, YSTR: Strings that are eval'd on the data returned from
%    quick_load to yield the x and y values of the data to be plotted.
%    For example, if XSTR were 'data(1,:,1)', then d.data(1,:,1) would be
%    used as the x-value of data to be plotted.
%
%    XFORM: A function that is applied to the data to be plotted prior to
%    mean and standard-error calculation.  This can be used to "unwind"
%    fitness functions, normalize, etc.  For example, @(i,x) (x./1600)
%    would divide all loaded y-values by 1600.
%
%    EXAMPLE:
%      d = quick_plot('data(1,:,1)', 'data(:,:,2)',...
%          {'001-structure/ini_.*gen-Genchamp-AvgFit', 'hn-ini';...
%           '001-structure/geo1_.*gen-Genchamp-AvgFit', 'hn-geo1';...
%           '001-structure/3d_.*gen-Genchamp-AvgFit', '3d';...
%           '001-structure/mhn_.*gen-Genchamp-AvgFit', 'ma-hn';...
%           '001-structure/n_.*gen-Genchamp-AvgFit', 'n'},...
%          @(i,x) (x./1600));
%      axis([0 200 0 1.05]);
%      xlabel('Generations')
%      ylabel('Fraction of max. fitness');
%      title('Varied neuroevolutionary approaches');
%      quick_export('../doc/gecco2010/figures/001-structure.eps');
if(nargin < 4); xform = @(i,x) x; end

markers_per_line=12;
marker_size=8;
line_width=1.0;
% color, marker, std err style
style = {'b','^',':';...
    'r','o',':';...
    'g','*',':';...
    'm','+',':';...
    'k','s',':';...
    'y','v',':'};


d = cell(size(ql,1),1);
for i=1:size(ql,1)
    d{i} = quick_load(ql{i,1});
end

% [lb0,ub0] = all_bootci(y0,200,@mean,0.05);
% [lb1,ub1] = all_bootci(y1,200,@mean,0.05);
% [lb2,ub2] = all_bootci(y2,200,@mean,0.05);
% [lb3,ub3] = all_bootci(y3,200,@mean,0.05);
% 
% [fig axes] = newfigure();
% ciplot(lb0(1:4:end), ub0(1:4:end), x(1:4:end), 'r');
% ciplot(lb1(1:4:end), ub1(1:4:end), x(1:4:end), 'b');
% ciplot(lb2(1:4:end), ub2(1:4:end), x(1:4:end), 'g');
% ciplot(lb3(1:4:end), ub3(1:4:end), x(1:4:end), 'y');
% plot(x, mean(y0), 'r', 'linewidth', 1.5);


[fig axes] = newfigure();
lh = [ ];
for i=1:size(d,1)
    x = eval(['d{' num2str(i) '}.' xstr]);
    y = eval(['d{' num2str(i) '}.' ystr]);
    y = xform(i,y);
    ym = mean(y);
    yse = stderr(y);
    skip = floor(size(x,2)/markers_per_line);
    
    p = plot(x, ym, style{i,1},...
        x(1:skip:end), ym(1:skip:end), [style{i,1} style{i,2}],...
        x, ym+yse, [style{i,1} style{i,3}],...
        x, ym-yse, [style{i,1} style{i,3}],...
        'markersize', marker_size, 'linewidth', line_width);
    lh(i) = p(1);
end

legend(lh, ql{:,2}, 'orientation','horizontal','location','southoutside');
[lh, oh, ph, s] = legend(axes);

for i=1:size(d,1)
    set(oh(size(d,1)+i*2),'marker',style{i,2});
end