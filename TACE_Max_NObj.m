function [maxima, thresholded] = TACE_Max_NObj(img, n)
    % TACE_Max_NObj - Finds n brightest objects in descending intensity order
    % Enhanced version with improved peak detection and sorting
    
    % Input Validation
    if nargin < 2
        error('Function requires both image and n as inputs');
    end
    
    if ~ismatrix(img)
        error('Input must be a 2D image matrix');
    end
    
    if n <= 0 || n ~= round(n)
        error('n must be a positive integer');
    end
    
    % Enhanced Preprocessing
    img = double(img);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)) + eps); % Add eps to avoid division by zero
    
    % Multi-scale filtering for better peak preservation
    sigma1 = 1.0;
    sigma2 = 2.5;
    img_filt1 = imgaussfilt(img, sigma1);
    img_filt2 = imgaussfilt(img, sigma2);
    filtered_img = 0.6*img_filt1 + 0.4*img_filt2;
    
    % Advanced background subtraction
    se = strel('disk', round(min(size(img))/20));
    background = imopen(filtered_img, se);
    subtracted_img = filtered_img - background;
    
    % Improved Peak Detection
    % Find all regional maxima
    conn = 8; % Connectivity for more precise peak detection
    regional_max = imregionalmax(subtracted_img, conn);
    
    % Get all peak positions and intensities
    [rows, cols] = find(regional_max);
    intensities = subtracted_img(regional_max);
    
    % Sort peaks by intensity (descending)
    [sorted_intensities, sort_idx] = sort(intensities, 'descend');
    sorted_positions = [rows(sort_idx), cols(sort_idx)];
    
    % Cluster Nearby Peaks (avoid duplicates)
    min_dist = 10; % Minimum distance between peaks
    maxima = [];
    current_idx = 1;
    
    while ~isempty(sorted_positions) && size(maxima,1) < n
        % Take brightest remaining peak
        current_peak = sorted_positions(1,:);
        maxima(current_idx,:) = current_peak;
        current_idx = current_idx + 1;
        
        % Remove peaks too close to this one
        distances = sqrt(sum((sorted_positions - current_peak).^2, 2));
        keep_idx = distances > min_dist;
        sorted_positions = sorted_positions(keep_idx,:);
        sorted_intensities = sorted_intensities(keep_idx);
    end
    
    % If not enough peaks found, use remaining sorted positions
    if size(maxima,1) < n
        needed = n - size(maxima,1);
        maxima = [maxima; sorted_positions(1:min(needed,end),:)];
    end
    
    % Enhanced Visualization
    thresholded = repmat(subtracted_img, [1 1 3]);
    thresholded = imadjust(thresholded, stretchlim(thresholded, [0.01 0.99]));
    
    % Mark maxima with clear indicators
    marker_size = max(5, ceil(min(size(img))/50));
    text_offset = marker_size * 2;
    
    for i = 1:size(maxima,1)
        % Draw numbered markers
        thresholded = insertShape(thresholded, 'FilledCircle', ...
            [maxima(i,2), maxima(i,1), marker_size], ...
            'Color', 'red', 'Opacity', 0.7);
        
        thresholded = insertText(thresholded, ...
            [maxima(i,2)+text_offset, maxima(i,1)-text_offset], ...
            i, 'TextColor', 'white', 'FontSize', 14, ...
            'BoxColor', 'red', 'BoxOpacity', 0.6);
    end
    
    % Optional Display
    if nargout == 0
        figure('Name','Optimized Bright Object Detection');
        subplot(1,3,1);
        imshow(img, []);
        title('Original Image');
        
        subplot(1,3,2);
        imshow(subtracted_img, []);
        title('Processed Image');
        
        subplot(1,3,3);
        imshow(thresholded);
        title(sprintf('Top %d Brightest Objects', n));
    end
end