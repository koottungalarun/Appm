function A = readSparseMatrix(filename)
% READSPARSEMATRIX    Read a sparse matrix from text file.
%
% The data file contains a sparse matrix in coordinate format, i.e., for
% each nonzero value, we have a triplet (row, column, value).
%
% Expected file format: (nnz rows, 3 columns)
% 
% r, c, v
% r, c, v
% (...)
%
% Matlab requires a one-based format.
% If a row/column index is equal to zero, then we assume a zero-based index
% format; consequently, this is converted to one-based format.
%
% Author: Roman Fuchs
% Date: 15.04.2019
%



%% Check input parameters
if isempty(filename)
    error('Input parameter for filename is empty.');
end

if ~exist(filename, 'file')
    error('File does not exist: %s', filename);
end

%% Read data file
data = importdata(filename);
if isempty(data)
    A = [];
    return
end
assert(size(data,2) == 3);
assert(size(data,1) > 0);

% Determine file format: zero-based or one-based?
isZeroBased = any((data(:,1) == 0)) || any(data(:,2) == 0);
isOneBased = all(data(:,1) >= 1) && all(data(:,2) >= 1);

% File format should be unique, it can't be both or none of them
assert(isZeroBased ~= isOneBased);

% If file format is zero-based, increase indices by 1. 
% This ensures that we have a one-based format
if isZeroBased
    data(:,1:2) = data(:,1:2) + 1;
end


%% Create sparse matrix from file contents
% Matlab requires matrices to be in one-based format
A = spconvert(data);


end