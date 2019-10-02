function d = quick_multiload(ql, varargin)

d = cell(size(ql,1),1);
for i=1:size(ql,1)
    d{i} = quick_load(ql{i,1}, varargin{:});
end