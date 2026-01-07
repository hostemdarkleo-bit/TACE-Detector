function [SMatrix, MtPxData] = SMosaics(B, s, t)
% SMOSAICS Divide an MxN matrix into 's' equal mosaics (tiles).
% 
% This function partitions large astronomical images into smaller tiles for
% parallel processing in the topological variational snake algorithm.
% The tiling approach enables efficient processing of high-resolution
% images by dividing them into manageable sub-regions.
%
% Reference: "Topological_Variational_Snake.pdf"
%
% Inputs:
%   B - Input matrix of size MxN (M, N > 10000 for astronomical images).
%   s - Total number of mosaics (integer, s >= 4).
%   t - Number of mosaic rows (integer, t >= 2 and must be a divisor of 's').
%
% Outputs:
%   SMatrix    - 3D matrix of size mxnxs containing all mosaics.
%   MtPxData   - Structure array containing position and size information for each mosaic.
%
% Usage in Topological Variational Snake Algorithm:
%   This function enables the division of large astronomical images into
%   smaller regions that can be processed independently, facilitating:
%   1. Parallel computation of object detection
%   2. Memory-efficient handling of high-resolution images
%   3. Independent optimization of variational contours in each region

% --- INPUT VALIDATION ---
[M, N] = size(B);

% Validate 's': must be an integer and greater than or equal to 4
% This ensures meaningful partitioning of the image
if ~isscalar(s) || ~isnumeric(s) || s < 4 || floor(s) ~= s
    error('SMosaics:InvalidS', ...
          'The value of ''s'' (total number of mosaics) must be a single integer greater than or equal to 4.');
end

% Validate 't': must be an integer, greater than or equal to 2, and a divisor of 's'
% This ensures proper rectangular tiling arrangement
if ~isscalar(t) || ~isnumeric(t) || t < 2 || floor(t) ~= t || mod(s, t) ~= 0
    error('SMosaics:InvalidT', ...
          'The value of ''t'' (number of mosaic rows) must be a single integer, greater than or equal to 2, and a divisor of ''s''.');
end

% --- DIMENSION CALCULATION AND PARTITIONING ---

% Determine mosaic arrangement based on t (rows) and s/t (columns)
% This creates a rectangular grid of tiles
rows = t;
cols = s / t;

% Calculate the dimensions of each mosaic (m x n)
% Using floor ensures integer dimensions; any remainder at the edges
% of the original image will be excluded from the mosaics
% For astronomical images, the excluded portion is typically negligible
m = floor(M / rows); % Height of each mosaic
n = floor(N / cols); % Width of each mosaic

% Check if the calculated dimensions are valid
% If m or n is less than 1, the image is too small for the requested partitioning
if m < 1 || n < 1
    error('SMosaics:MatrixTooSmall', ...
          'The dimensions of the input matrix B are too small for the requested number of mosaics (s and t).');
end

% Initialize the 3D matrix to store all mosaics
% Pre-allocating memory is crucial for efficiency with large astronomical images
% Using the same class as B preserves the original data type (e.g., double, uint16)
SMatrix = zeros(m, n, s, class(B));

% Initialize the structure array for mosaic metadata
% Pre-allocating the structure improves performance
MtPxData(s).position = []; % Pre-allocate to final size by setting last element
MtPxData(s).size = [];     % This establishes the structure array dimensions

% --- EXTRACTION AND STORAGE OF MOSAICS ---
tileCount = 0; % Use zero-based counter for loop indices

% Iterate through rows and columns to extract each mosaic
% This nested loop generates mosaics in row-major order (left to right, top to bottom)
for i = 0:(rows - 1) % Iterate through mosaic rows
    for j = 0:(cols - 1) % Iterate through mosaic columns
        tileCount = tileCount + 1; % Increment counter for current mosaic index
        
        % Calculate row and column boundaries for the current mosaic
        % These indices define the region of the original image to extract
        rowStart = i * m + 1;
        rowEnd = (i + 1) * m;
        colStart = j * n + 1;
        colEnd = (j + 1) * n;
        
        % Extract the mosaic from the original matrix B
        % Store it in the corresponding layer of the 3D SMatrix
        SMatrix(:,:,tileCount) = B(rowStart:rowEnd, colStart:colEnd);
        
        % Store metadata for the current mosaic
        % Position: [start_row, end_row, start_column, end_column]
        % Size: [height, width] of the mosaic
        % This information is crucial for mapping results back to the original image coordinates
        MtPxData(tileCount).position = [rowStart, rowEnd, colStart, colEnd];
        MtPxData(tileCount).size = [m, n];
    end
end

% End of function
end