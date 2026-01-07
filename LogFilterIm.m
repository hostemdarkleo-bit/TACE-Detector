function LogFiltImage = LogFilterIm(inputMatrix)
% LOGFILTERIM Applies an interactive logarithmic filter to an MxN matrix.
% 
% This function provides interactive adjustment of logarithmic contrast
% enhancement through a slider interface, allowing real-time visualization
% and saving of filtered images. The logarithmic transformation is
% particularly effective for enhancing low-intensity features in astronomical
% images while compressing high-intensity regions.
%
% Mathematical operation: g(x,y) = C * log(1 + f(x,y))
% where f(x,y) is the input image in [0,1] and C is a scaling constant
%
% Reference: "Topological_Variational_Snake.pdf"
%
% Inputs:
%   inputMatrix - The mxn sub-matrix to which the filter will be applied (double type).
%                 Typically this is a mosaic or sub-region of an astronomical image.
%
% Outputs:
%   LogFiltImage - The last filtered image when the window is closed.
%                  This will be the image displayed at the moment of closing the figure,
%                  or the last saved image if the "Save" button was used.
%
% Usage in Topological Variational Snake Algorithm:
%   This function enhances contrast in astronomical images to improve
%   the detection of faint objects and features that may be obscured
%   by background noise or uneven illumination.

% --- 1. PREPARE INPUT IMAGE ---
% Convert input matrix to double precision in range [0, 1]
% This normalization is essential for the logarithmic transformation
originalImage = im2double(inputMatrix);

% Variable to store the last filtered image (initialize with original)
LogFiltImage = originalImage;

% --- 2. CREATE FIGURE AND AXES FOR IMAGE DISPLAY ---
hFig = figure('Name', 'Logarithmic Filter', 'NumberTitle', 'off', 'Position', [100 100 800 650]);
hAx = axes('Parent', hFig, 'Position', [0.05 0.25 0.9 0.65]);

% Display the original image initially
hImg = imshow(originalImage, 'Parent', hAx);
title(hAx, 'Original Image');
axis image off;

% --- 3. CREATE LOGARITHMIC SCALING CONSTANT SLIDER ---
% Define slider parameters for logarithmic scaling constant C
% The constant C controls the strength of the logarithmic transformation
minC = 1;          % Minimum scaling constant (weak enhancement)
maxC = 100;        % Maximum scaling constant (strong enhancement)
initialC = 10;     % Initial value (moderate enhancement)

% Create slider control for interactive adjustment of logarithmic constant
hSlider = uicontrol('Parent', hFig, 'Style', 'slider', ...
                    'Min', minC, 'Max', maxC, ...
                    'Value', initialC, ...
                    'Position', [150 100 500 20], ...
                    'Callback', @sliderCallback);

% Create text label to display current logarithmic constant
hText = uicontrol('Parent', hFig, 'Style', 'text', ...
                  'Position', [350 125 100 20], ...
                  'String', sprintf('Constant C: %.1f', initialC));

% --- 4. CREATE SAVE BUTTON ---
% Button to save the current filtered image to file
hSaveButton = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
                        'String', 'Save Image', ...
                        'Position', [350 30 100 30], ...
                        'Callback', @saveButtonCallback);

% --- 5. DEFINE CALLBACK FUNCTIONS ---

    % Slider callback: updates image display when slider value changes
    function sliderCallback(~, ~)
        % Get current slider value (logarithmic scaling constant C)
        currentC = get(hSlider, 'Value');
        
        % Apply logarithmic transformation: g(x,y) = C * log(1 + f(x,y))
        % The input image is in range [0,1], making log(1+x) appropriate
        % This transformation enhances low-intensity features while
        % compressing high-intensity regions
        tempFilteredImage = currentC * log(1 + originalImage);
        
        % Normalize the resulting image for display
        % 'imadjust' with 'stretchlim' maps values to [0, 1] range
        % This ensures optimal contrast for visualization
        if max(tempFilteredImage(:)) > 0
            LogFiltImage = imadjust(tempFilteredImage, stretchlim(tempFilteredImage), []);
        else 
            % Handle edge case: if image has no variation or is entirely dark
            LogFiltImage = tempFilteredImage;
        end

        % Update displayed image with filtered result
        set(hImg, 'CData', LogFiltImage);
        
        % Update slider text label
        set(hText, 'String', sprintf('Constant C: %.1f', currentC));
        
        % Update plot title
        title(hAx, sprintf('Logarithmic Filter (C: %.1f)', currentC));
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
            % Note: LogFiltImage is already normalized to [0,1] by imadjust
            imwrite(LogFiltImage, fullFileName);
            
            % Display confirmation message in command window
            disp(['Image saved as: ', fullFileName]);
        end
    end

% --- 6. FIGURE CLOSE CALLBACK ---
% Set callback for when figure is closed to ensure output variable contains last image
set(hFig, 'CloseRequestFcn', @figureCloseCallback);

    function figureCloseCallback(~,~)
        % The 'LogFiltImage' variable already contains the last displayed image
        % thanks to its update in sliderCallback
        delete(hFig); % Close the figure
    end

% --- 7. INITIAL CALL AND WAIT ---
% Call slider callback initially to display the image with default constant
sliderCallback([], []);

% Optional: Wait for the figure to close before returning the value
% This ensures the function returns the final filtered image
uiwait(hFig);

% End of function
end