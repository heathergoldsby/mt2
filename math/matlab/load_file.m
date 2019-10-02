function d = load_file(file, format)
%LOAD_FILE    Loads data from a single file.
%
%   d = LOAD_FILE(file, format) loads data from the specified file
%   according to the given format.  Automatically handles gzipped files.
%   See textscan for a description of the format field.
%
%   The returned value, d, is an m x n matrix for the rows and columns in
%   the file, respectively.

[status result] = unix(['file ' file]);
if status ~= 0
    error(['(load_file) error determining file type while loading' file]);
end

if regexp(result, 'gzip')
    [status result] = unix(['gzcat ' file]);
    if status ~= 0
        error(['(load_file) error in gzcat while loading' file]);
    end
else
    [status result] = unix(['cat ' file]);
    if status ~= 0
        error(['(load_file) error in cat while loading' file]);
    end
end

d = textscan(result, format, 'CommentStyle', '#');
