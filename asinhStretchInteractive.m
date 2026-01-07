function asinhStretchInteractive(inputArg)
% ASINHSTRETCHINTERACTIVE
% Interactive asinh stretch with sliders.
% Supports:
% - Workspace variable (matrix or variable name)
% - Files: png, jpg, tif, gif, fit, fits, m
%
% Usage:
%   asinhStretchInteractive()
%   asinhStretchInteractive(I)
%   asinhStretchInteractive("I")
%   asinhStretchInteractive("C:\path\image.fits")

% This function is part of the semi-automatic model based on topological
% isophotes in Sobolev spaces for object identification in astrophysical images.
% It provides an interactive interface to apply asinh stretching to images,
% which enhances low-intensity features while preserving high-intensity details.
% Reference: "Topological_Variational_Snake.pdf"

    % If no input is provided, prompt user to select input source
    if nargin < 1
        choice = questdlg('Select input source:', 'Asinh Stretch', ...
                          'Workspace variable','File','Workspace variable');
        if isempty(choice), return; end
        if strcmp(choice,'Workspace variable')
            I = pickWorkspaceImage();
        else
            I = pickFileImage();
        end
    else
        % Handle different input types
        if isnumeric(inputArg) || islogical(inputArg)
            I = toRGB01(inputArg);
        elseif isstring(inputArg) || ischar(inputArg)
            s = char(inputArg);
            if isfile(s)
                I = readInputFlexible(s);
            else
                I = toRGB01(evalin('base', s));
            end
        else
            error('Unsupported input type.');
        end
    end

    if isempty(I), return; end

    % ---- UI CREATION ----
    % Create main figure for the interactive interface
    fig = uifigure('Name','Interactive Asinh Stretch (R2024b)', ...
                   'Position',[100 100 1200 720]);

    % Create grid layout for image display and control panel
    gl = uigridlayout(fig,[1 2]);
    gl.ColumnWidth = {'1x', 360};
    gl.RowHeight = {'1x'};

    % Create axes for image display
    ax = uiaxes(gl);
    ax.Toolbar.Visible = 'off';
    ax.XTick = []; ax.YTick = [];
    ax.Box = 'on';

    % Create control panel for parameters
    pRight = uipanel(gl,'Title','Parameters','FontWeight','bold');
    gr = uigridlayout(pRight,[14 2]);
    gr.RowHeight = {22, 28, 22, 28, 22, 28, 22, 28, 22, 28, 22, 28, 28, '1x'};
    gr.ColumnWidth = {140,'1x'};

    % Create UI controls for stretch parameters
    uilabel(gr,'Text','Mode','FontWeight','bold');
    ddMode = uidropdown(gr,'Items',{'Luminance (recommended)','Per-channel (RGB)'}, ...
                           'Value','Luminance (recommended)');

    uilabel(gr,'Text','Black point B');
    sB = uislider(gr,'Limits',[0 1],'Value',0.02,'MajorTicks',[0 0.25 0.5 0.75 1]);

    uilabel(gr,'Text','White point W');
    sW = uislider(gr,'Limits',[0 1],'Value',0.98,'MajorTicks',[0 0.25 0.5 0.75 1]);

    uilabel(gr,'Text','Strength beta');
    sBeta = uislider(gr,'Limits',[0.1 500],'Value',20,'MajorTicks',[0.1 10 50 200 500]);

    uilabel(gr,'Text','Offset eps');
    sE = uislider(gr,'Limits',[0 0.05],'Value',0.002,'MajorTicks',[0 0.01 0.02 0.03 0.04 0.05]);

    uilabel(gr,'Text','Gamma (output)');
    sG = uislider(gr,'Limits',[0.3 3],'Value',1.0,'MajorTicks',[0.3 0.7 1 1.5 2 3]);

    uilabel(gr,'Text','Options','FontWeight','bold');
    cbHist = uicheckbox(gr,'Text','Show luminance histogram','Value',false);

    % Create control buttons
    btnReset = uibutton(gr,'Text','Reset','ButtonPushedFcn',@(~,~)resetParams()); %#ok<NASGU>
    btnSave  = uibutton(gr,'Text','Save PNG','ButtonPushedFcn',@(~,~)saveResult()); %#ok<NASGU>

    % Create label for displaying current parameter values
    lblVals = uilabel(gr,'Text','','VerticalAlignment','top');
    lblVals.Layout.Row = [13 14];
    lblVals.Layout.Column = [1 2];

    % Initialize variables for histogram display
    hFigHist = [];
    hAxHist  = [];
    hLineHist = [];

    % Display the original image
    hImg = imshow(I,'Parent',ax);
    title(ax,'Loaded input','FontWeight','bold');
    updateDisplay();

    % Set up callback functions for interactive parameter adjustment
    % ValueChangingFcn provides real-time updates during slider movement
    sB.ValueChangingFcn = @(~,evt) updateDisplay(evt.Value, []);
    sW.ValueChangingFcn = @(~,evt) updateDisplay([], evt.Value);

    % ValueChangedFcn provides updates after slider release
    sBeta.ValueChangingFcn = @(~,~) updateDisplay();
    sE.ValueChangingFcn    = @(~,~) updateDisplay();
    sG.ValueChangingFcn    = @(~,~) updateDisplay();

    sB.ValueChangedFcn = @(~,~) updateDisplay();
    sW.ValueChangedFcn = @(~,~) updateDisplay();
    sBeta.ValueChangedFcn = @(~,~) updateDisplay();
    sE.ValueChangedFcn = @(~,~) updateDisplay();
    sG.ValueChangedFcn = @(~,~) updateDisplay();

    ddMode.ValueChangedFcn = @(~,~) updateDisplay();
    cbHist.ValueChangedFcn = @(~,~) updateDisplay();

    % ---------------- NESTED FUNCTIONS ----------------

    % Main display update function - called whenever parameters change
    function updateDisplay(Btemp, Wtemp)
        % Get current parameter values, allowing for temporary values during slider drag
        B = sB.Value;
        W = sW.Value;
        if nargin >= 1 && ~isempty(Btemp), B = Btemp; end
        if nargin >= 2 && ~isempty(Wtemp), W = Wtemp; end
        
        % Ensure white point is greater than black point
        if W <= B + 1e-6
            W = min(1, B + 1e-3);
        end

        beta = sBeta.Value;
        eps0 = sE.Value;
        gammaOut = sG.Value;
        mode = ddMode.Value;

        % Apply asinh stretch and gamma correction
        Iout = applyAsinhStretch(I, B, W, beta, eps0, mode);
        Iout = Iout .^ (1/gammaOut);

        % Update displayed image
        hImg.CData = Iout;
        title(ax, sprintf('Asinh stretch | B=%.4f W=%.4f beta=%.2f eps=%.4f gamma=%.2f', ...
            B, W, beta, eps0, gammaOut), 'FontWeight','bold');

        % Update parameter display label
        lblVals.Text = sprintf([ ...
            'Tip: tune B/W first, then beta, then eps, then gamma.\n' ...
            'Mode: %s\n' ...
            'B=%.4f | W=%.4f | beta=%.2f | eps=%.4f | gamma=%.2f'], ...
            mode, B, W, beta, eps0, gammaOut);

        % Show/hide histogram based on checkbox state
        if cbHist.Value
            L = getLuminance(Iout);
            showHistogram(L);
        else
            if ~isempty(hFigHist) && isvalid(hFigHist), close(hFigHist); end
            hFigHist = [];
        end

        drawnow limitrate;
    end

    % Core asinh stretching function
    function Iout = applyAsinhStretch(Iin, B, W, beta, eps0, mode)
        % Normalize input to range [0,1] based on black and white points
        X = (Iin - B) / (W - B);
        X = min(max(X,0),1);

        % Pre-compute denominator for asinh normalization
        denom = asinh(beta*(1 + eps0));

        if contains(mode,'Luminance')
            % Luminance-based stretching: apply stretch to luminance channel only
            L  = getLuminance(X);
            Ls = asinh(beta*(L + eps0)) / denom;
            Ls = min(max(Ls,0),1);

            % Scale RGB channels by luminance ratio
            tiny = 1e-6;
            scale = Ls ./ (L + tiny);
            Iout = X .* scale;
            Iout = min(max(Iout,0),1);
        else
            % Per-channel stretching: apply stretch independently to each RGB channel
            Iout = asinh(beta*(X + eps0)) / denom;
            Iout = min(max(Iout,0),1);
        end
    end

    % Calculate luminance from RGB image using standard coefficients
    function L = getLuminance(Irgb)
        L = 0.2126*Irgb(:,:,1) + 0.7152*Irgb(:,:,2) + 0.0722*Irgb(:,:,3);
    end

    % Display histogram of luminance values
    function showHistogram(L)
        % Create histogram figure if it doesn't exist
        if isempty(hFigHist) || ~isvalid(hFigHist)
            hFigHist = uifigure('Name','Luminance Histogram','Position',[1350 200 520 360]);
            hAxHist = uiaxes(hFigHist);
            hAxHist.Box = 'on';
            hAxHist.XLim = [0 1];
            title(hAxHist,'Output luminance histogram');
            xlabel(hAxHist,'Intensity');
            ylabel(hAxHist,'Count');
        end
        
        % Compute and display histogram
        nbins = 256;
        [counts, edges] = histcounts(L(:), nbins, 'BinLimits',[0 1]);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        
        % Update or create histogram plot
        if isempty(hLineHist) || ~isvalid(hLineHist)
            hLineHist = plot(hAxHist, centers, counts, 'LineWidth', 1.2);
        else
            hLineHist.XData = centers;
            hLineHist.YData = counts;
        end
    end

    % Reset all parameters to default values
    function resetParams()
        sB.Value = 0.02;
        sW.Value = 0.98;
        sBeta.Value = 20;
        sE.Value = 0.002;
        sG.Value = 1.0;
        ddMode.Value = 'Luminance (recommended)';
        cbHist.Value = false;
        updateDisplay();
    end

    % Save the processed image to file
    function saveResult()
        % Get current parameter values
        B = sB.Value; W = sW.Value;
        if W <= B + 1e-6, W = min(1, B + 1e-3); end
        beta = sBeta.Value;
        eps0 = sE.Value;
        gammaOut = sG.Value;
        mode = ddMode.Value;

        % Apply stretch with current parameters
        Iout = applyAsinhStretch(I, B, W, beta, eps0, mode);
        Iout = Iout .^ (1/gammaOut);

        % Prompt user for save location
        [f2,p2] = uiputfile('asinhstretch_output.png','Save result as');
        if isequal(f2,0), return; end
        
        % Save image and notify user
        imwrite(Iout, fullfile(p2,f2));
        uialert(fig,'Image saved successfully.','OK');

        % Export to workspace for use in subsequent processing
        assignin('base','LogFilt_image', Iout);
    end

end

% ---------------- HELPER FUNCTIONS ----------------

% Select image file via file dialog
function I = pickFileImage()
    [f,p] = uigetfile({ ...
        '*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.gif;*.fit;*.fits;*.m', ...
        'Supported files (*.png, *.jpg, *.tif, *.gif, *.fit, *.fits, *.m)'; ...
        '*.*','All files (*.*)'}, ...
        'Select file');
    if isequal(f,0), I = []; return; end
    I = readInputFlexible(fullfile(p,f));
end

% Select image from workspace variables
function I = pickWorkspaceImage()
    % Get list of variables from base workspace
    ws = evalin('base','whos');
    names = {};
    for k = 1:numel(ws)
        % Filter for numeric/image-like variables
        if (strcmp(ws(k).class,'double') || strcmp(ws(k).class,'single') || ...
            strcmp(ws(k).class,'uint8')  || strcmp(ws(k).class,'uint16') || ...
            strcmp(ws(k).class,'int16')  || strcmp(ws(k).class,'logical')) ...
            && numel(ws(k).size) >= 2
            names{end+1} = ws(k).name; %#ok<AGROW>
        end
    end
    
    if isempty(names)
        uialert(uifigure,'No numeric image-like variables found in base workspace.','Workspace');
        I = [];
        return;
    end
    
    % Let user select from available variables
    [idx, ok] = listdlg('PromptString','Select a variable from base workspace:', ...
                        'SelectionMode','single', 'ListString', names);
    if ~ok, I = []; return; end
    
    A = evalin('base', names{idx});
    I = toRGB01(A);
end

% Flexible image reader supporting multiple file formats
function I = readInputFlexible(filePath)
    [~,~,ext] = fileparts(filePath);
    ext = lower(ext);

    % Handle FITS files (common in astronomy)
    if any(strcmp(ext,{'.fit','.fits'}))
        A = readFITSFirstPlane(filePath);
        I = toRGB01(A);
        return;
    end

    % Handle GIF files with color maps
    if strcmp(ext,'.gif')
        [A, map] = imread(filePath, 1);
        if ~isempty(map)
            A = ind2rgb(A, map);
            I = im2double(A);
        else
            I = im2double(A);
        end
        I = ensureRGB(I);
        return;
    end

    % Handle MATLAB .m files that may contain image data
    if strcmp(ext,'.m')
        warning('Executing M-file to load image data. Use only trusted .m files.');
        A = readFromMFileTrusted(filePath);
        I = toRGB01(A);
        return;
    end

    % Handle standard image formats
    A = imread(filePath);
    I = im2double(A);
    I = ensureRGB(I);
end

% Read first plane of FITS file (astronomical data format)
function A = readFITSFirstPlane(filePath)
    try
        try
            A = fitsread(filePath,'primary');
        catch
            A = fitsread(filePath);
        end
    catch ME
        error('Could not read FITS: %s', ME.message);
    end
    
    A = double(A);
    A(~isfinite(A)) = 0;

    % Handle multi-dimensional FITS data
    sz = size(A);
    if numel(sz) > 3
        idx = repmat({':'},1,numel(sz));
        for k = 4:numel(sz), idx{k} = 1; end
        A = A(idx{:});
    end
    
    % Extract first channel if not RGB
    if ndims(A) == 3 && size(A,3) ~= 3
        A = A(:,:,1);
    end
end

% Read image data from MATLAB .m file
function A = readFromMFileTrusted(filePath)
    S = runMfileToStruct(filePath);
    
    % Try common variable names for image data
    candidates = {'I','img','image','A','data','frame'};
    for k = 1:numel(candidates)
        nm = candidates{k};
        if isfield(S,nm) && (isnumeric(S.(nm)) || islogical(S.(nm)))
            A = S.(nm);
            return;
        end
    end
    
    % If no common names found, let user select from available variables
    fn = fieldnames(S);
    keep = false(size(fn));
    for i = 1:numel(fn)
        v = S.(fn{i});
        keep(i) = isnumeric(v) || islogical(v);
    end
    fn = fn(keep);
    
    if isempty(fn), error('No numeric variables found in the M-file.'); end
    
    [idx, ok] = listdlg('PromptString','Select the variable that contains the image matrix:', ...
                        'SelectionMode','single', 'ListString', fn);
    if ~ok, error('No variable selected.'); end
    
    A = S.(fn{idx});
end

% Execute .m file and return all variables as struct
function S = runMfileToStruct(filePath)
    function out = runner()
        run(filePath);
        w = whos;
        out = struct();
        for i = 1:numel(w)
            name = w(i).name;
            out.(name) = eval(name);
        end
    end
    S = runner();
end

% Convert any numeric array to normalized RGB image in [0,1] range
function I = toRGB01(A)
    A = double(A);
    A(~isfinite(A)) = 0;
    
    % If already RGB, normalize each channel
    if ndims(A) == 3 && size(A,3) == 3
        I = zeros(size(A));
        for c = 1:3
            I(:,:,c) = normalizeTo01(A(:,:,c));
        end
    else
        % If grayscale, convert to RGB
        if ndims(A) == 3, A = A(:,:,1); end
        G = normalizeTo01(A);
        I = repmat(G,1,1,3);
    end
    I = min(max(I,0),1);
end

% Normalize array to [0,1] range using percentile-based scaling
function X = normalizeTo01(A)
    v = A(:);
    v = v(isfinite(v));
    if isempty(v), X = zeros(size(A)); return; end
    
    % Check if data is already in [0,1] range
    if prctile(v,99.9) <= 1.2 && prctile(v,0.1) >= -0.2
        X = min(max(A,0),1);
        return;
    end
    
    % Use percentile-based normalization for robustness to outliers
    pLow  = prctile(v, 0.5);
    pHigh = prctile(v, 99.5);
    if pHigh <= pLow
        pLow = min(v); pHigh = max(v);
        if pHigh <= pLow, X = zeros(size(A)); return; end
    end
    
    X = (A - pLow) / (pHigh - pLow);
    X = min(max(X,0),1);
end

% Ensure image has 3 channels (RGB)
function I = ensureRGB(I)
    if ndims(I) == 2
        I = repmat(I,1,1,3);
    end
end