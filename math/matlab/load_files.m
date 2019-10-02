function [data,fnames] = load_files(filenames, format, keepall)
%LOAD_FILES    Loads multiple data files into a single matrix.

if(nargin < 2)
    error('(load_files) wrong number of input arguments');
end

% load the data from the files into a cell array:
x=cell(1,size(filenames,2));
for i=1:size(filenames,2), j=char(filenames{i});
    x{i} = load_file(j, format);
end

% calculate the ultimate dimensions of the returned struct
%
m=size(x,2); % each row belongs to a different trial
n=0; % each column is a single time
p=size(x{1},2); % each page is a different variable

% find the maximum number of times in our data:
for i=1:m;
    if size(x{i}{1},1) > n
        n = size(x{i}{1},1);
    end
end

% now, which of our trials have this data?
trials=ones(m,1); % logical
ctrials=m;
if ~keepall
    for i=1:m;
        if size(x{i}{1},1) < n
            disp(['Warning: file ' char(filenames{i}) ' appears to be missing data.']);
            trials(i) = false;
            ctrials = ctrials - 1;
        end
    end
end

% preallocate data
data = NaN(ctrials,n,p);
fnames = cell(1,ctrials);

% and now copy the data and associated filenames that we're keeping:
cdata=1;
for i=1:m
    if trials(i)
        fnames{cdata} = filenames{i};
        for j=1:p;
            data(cdata, 1:size(x{i}{1},1), j) = x{i}{j};
        end
        cdata = cdata + 1;
    end
end
