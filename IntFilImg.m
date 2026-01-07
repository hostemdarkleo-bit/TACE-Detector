function AtImage = IntFilImg(inputMatrix)
% INTFILIMG Applies an interactive intensity reduction filter to an MxN matrix.
% 
% This function provides interactive adjustment of image intensity through
% a slider interface, allowing real-time visualization and saving of filtered
% images. It is part of the preprocessing pipeline for the topological
% variational snake algorithm used in astrophysical object detection.
%
% Reference: "Topological_Variational_Snake.pdf"
%
% Inputs:
%   inputMatrix - The mxn sub-matrix to which the filter will be applied (double type).
%                 Typically this is a mosaic or sub-region of an astronomical image.
%
% Outputs:
%   AtImage - The last filtered image when the window is closed.
%             This will be the image displayed at the moment of closing the figure,
%             or the last saved image if the "Save" button was used.
%
% Usage in Topological Variational Snake Algorithm:
%   This function enables manual adjustment of image intensity to optimize
%   object detection in challenging astronomical images where automatic
%   methods may fail due to varying background levels or noise.

% --- 1. PREPARE INPUT IMAGE ---
% Ensure the image is in double precision and normalized to [0, 1] range
% This is crucial for consistent handling if inputMatrix comes from a different range or data type
originalImage = im2double(inputMatrix);

% Variable to store the last filtered image (initialize with original)
AtImage = originalImage;

% --- 2. CREATE FIGURE AND AXES FOR IMAGE DISPLAY ---
hFig = figure('Name', 'Intensity Reduction Filter', 'NumberTitle', 'off', 'Position', [100 100 800 650]);
hAx = axes('Parent', hFig, 'Position', [0.05 0.25 0.9 0.65]); % Adjust position for sliders and buttons

% Display the original image initially
hImg = imshow(originalImage, 'Parent', hAx);
title(hAx, 'Original Image');
axis image off;

% --- 3. CREATE INTENSITY SLIDER ---
% Define slider parameters for intensity adjustment
% The intensity factor scales pixel values: output = input × factor
minIntensityFactor = 0.1;   % Minimum scaling factor (10% intensity)
maxIntensityFactor = 1.0;   % Maximum scaling factor (100% intensity - no reduction)
initialIntensityFactor = 1.0; % Initial value (no intensity reduction)

% Create slider control for interactive intensity adjustment
hSlider = uicontrol('Parent', hFig, 'Style', 'slider', ...
                    'Min', minIntensityFactor, 'Max', maxIntensityFactor, ...
                    'Value', initialIntensityFactor, ...
                    'Position', [150 100 500 20], ...
                    'Callback', @sliderCallback);

% Create text label to display current intensity factor
hText = uicontrol('Parent', hFig, 'Style', 'text', ...
                  'Position', [350 125 100 20], ...
                  'String', sprintf('Factor: %.2f', initialIntensityFactor));

% --- 4. CREATE SAVE BUTTON ---
% Button to save the current filtered image to file
hSaveButton = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
                        'String', 'Save Image', ...
                        'Position', [350 30 100 30], ...
                        'Callback', @saveButtonCallback);

% --- 5. DEFINE CALLBACK FUNCTIONS ---

    % Slider callback: updates image display when slider value changes
    function sliderCallback(~, ~)
        % Get current slider value (intensity factor)
        currentFactor = get(hSlider, 'Value');
        
        % Apply intensity reduction: scale original image by current factor
        AtImage = originalImage * currentFactor; % Update output variable
        
        % Update displayed image
        set(hImg, 'CData', AtImage);
        
        % Update slider text label
        set(hText, 'String', sprintf('Factor: %.2f', currentFactor));
        
        % Update plot title
        title(hAx, sprintf('Intensity Reduced (Factor: %.2f)', currentFactor));
    end

    % Save button callback: saves current filtered image to file
    function saveButtonCallback(~, ~)
        % Prompt user for save location and filename
        [filename, pathname] = uiputfile({'*.png';'*.jpg';'*.tif'}, 'Save Filtered Image');
        
        % Check if user didn't cancel the dialog
        if ischar(filename)
            % Construct full file path
            fullFileName = fullfile(pathname, filename);
            
            % Save the current filtered image
            % Note: imwrite automatically handles data type conversion
            imwrite(AtImage, fullFileName);
            
            % Display confirmation message in command window
            disp(['Image saved as: ', fullFileName]);
        end
    end

% --- 6. FIGURE CLOSE CALLBACK ---
% Set callback for when figure is closed to ensure output variable contains last image
set(hFig, 'CloseRequestFcn', @figureCloseCallback);

    function figureCloseCallback(~,~)
        % The 'AtImage' variable already contains the last displayed image
        % thanks to its update in sliderCallback
        delete(hFig); % Close the figure
    end

% --- 7. INITIAL CALL AND WAIT ---
% Call slider callback initially to display the image with default factor
sliderCallback([], []);

% Optional: Wait for the figure to close before returning the value
% This ensures the function returns the final filtered image
uiwait(hFig);

% End of function
end