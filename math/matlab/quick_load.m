function R = quick_load(file_pattern, varargin)
%QUICK_LOAD
%
dir='./';
format='';
keepall=0;
skip=1;

i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i},
        case 'dir'
            dir = varargin{i+1};
            i=i+1;
        case 'format'
            format = varargin{i+1};
            i=i+1;
        case 'keepall'
            keepall = varargin{i+1};
            i=i+1;
        otherwise
            argok=0;
    end
  end
  if ~argok, 
    disp(['(quick_load) bad argument: ' varargin{i}]);
  end
  i=i+1; 
end

% find the files we're loading, grab the header lines from one.
filenames = find_files(dir, file_pattern);
lines = header_lines(filenames);
lines = lines(~cellfun('isempty',regexp(lines,'^#\s+\d+[:\.]\s')));

% build the format string if it hasn't been specified
if isempty(format)
    format = arrayfun(@(x){'%n'},lines);
    format = cat(2,format{:});
end

% load the data
[ldata fnames] = load_files(filenames, format, keepall);

% build the returned struct.
tokens=regexp(lines,'.*\[([\w_]+)\]$','tokens');
R.filenames=fnames;
R.fieldnames=lines;
for i=1:size(ldata,3);
    R.(char(tokens{i}{:})) = ldata(:,:,i);
end    
