function [BinaryMask, IndividualMasks] = TACE_SegmentNObj(img, n, maxima)
    % TACE_SegmentNObj - Robust segmentation using maxima positions
    % Inputs:
    %   img    - Input image (2D matrix)
    %   n      - Number of objects to segment
    %   maxima - [n x 2] matrix of [row,col] positions
    %
    % Outputs:
    %   BinaryMask     - Combined binary mask of all objects
    %   IndividualMasks - Cell array of individual object masks

    %% Initialize
    img = double(img);
    img_norm = mat2gray(img); % Normalize to [0,1]
    BinaryMask = false(size(img));
    IndividualMasks = cell(n,1);
    
    %% Enhanced adaptive thresholding
    % Multi-scale smoothing for better thresholding
    img_smooth = imgaussfilt(img_norm, 2);
    
    % Global threshold with local adjustment
    global_thresh = graythresh(img_smooth);
    adaptive_thresh = imadjust(img_smooth, stretchlim(img_smooth), []);
    bw = imbinarize(adaptive_thresh, global_thresh*0.9); % More inclusive
    
    %% Process each maxima point
    for i = 1:min(n, size(maxima,1))
        try
            % Create ROI around maxima
            [rows, cols] = size(img);
            [xx,yy] = meshgrid(1:cols, 1:rows);
            
            % Adaptive radius estimation
            radius = estimate_radius(img_norm, maxima(i,:));
            
            % Create mask based on Snake. The Snake Level Set surfaces are
            % circle o elliptical curves.
            circle_mask = sqrt((xx-maxima(i,2)).^2 + (yy-maxima(i,1)).^2) <= radius;
            
            % Combine with thresholded image
            objMask = bw & circle_mask;
            
            % Get largest connected component
            objMask = bwareafilt(objMask, 1);
            objMask = imfill(objMask, 'holes');
            
            % Store results
            IndividualMasks{i} = objMask;
            BinaryMask = BinaryMask | objMask;
            
        catch ME
            warning('Error processing object %d: %s', i, ME.message);
            IndividualMasks{i} = false(size(img));
        end
    end
    
    %% Final cleanup
    BinaryMask = imfill(BinaryMask, 'holes');
    BinaryMask = bwareaopen(BinaryMask, 50); % Remove small artifacts
end
