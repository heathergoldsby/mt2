function quick_export(filename,w,h)
%QUICK_EXPORT    Export the current figure (gcf) to the specified file.
%
%    QUICK_EXPORT(filename,h,w) exports the current figure (gcf) to the
%    specified file in preferred publication-ready format.
if nargin < 3; w=6; h=3.75; end

exportfig(gcf, filename,...
        'Color','cmyk',...
        'Width', w, 'Height', h,...
        'FontMode', 'fixed', 'FontSize', 12);
