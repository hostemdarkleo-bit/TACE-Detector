%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code is a simple modification of DetTACE.
%    The modification consists of applying the "for" cycle to all
%    "s" mosaics from the original image.
%
%    This script extends the topological variational snake algorithm
%    to process multiple image mosaics simultaneously, applying object
%    detection and segmentation to each mosaic independently.
%    Reference: "Topological_Variational_Snake.pdf"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read metadata
%Img = fitsread('C:\MathematicsNew\Bolivia.fits');

% Load pre-processed image data from MATLAB file
% This approach is used when working with pre-computed image matrices
% rather than direct FITS file reading
image = load('C:\MathematicsNew\imagenes_el_gordo_Hubble\mats_el_gordo_hst\hst_mos_1041271_acs_wfc_f850lp_drz.mat');
Image = image.data;

% Launch interactive asinh stretching tool for image enhancement
% This allows dynamic adjustment of contrast and brightness for optimal
% object detection in astronomical images
asinhStretchInteractive()

% Convert the filtered RGB image to grayscale for processing
% The asinh stretching function returns an RGB image, but subsequent
% topological analysis requires grayscale input
Img = rgb2gray(LogFilt_image);

% Display the original study image
figure;
imshow(Img)
title('Original study Image');

% Image attenuation based only on intensity with attenuation factor "AtFac = 12"
% This step normalizes the image intensity range for consistent processing
% Note: AtFac=1 means no attenuation is applied
AtFac = 1;
Img = Img/AtFac;

% Partition into s tiles, arranged in t rows
% The original image is divided into s sub-images (mosaics)
% arranged in t rows for systematic processing
% This tiling approach allows parallel processing of different regions
% and handling of large astronomical images
s = 9;
t = 3;

% Number of objects of interest based on brightness per sub-image in the mosaic
% This parameter controls how many of the brightest objects will be
% detected and segmented in each mosaic
n = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partition where mosaics are aligned and arranged in descending order
% from the first entry of the top row to the last row from beginning to end
[SMatrix1, MtPxData1] = SMosaics(Img, s, t);

% Visualize tile positions
% Display all mosaics in a single figure for quality control
figure;
for i = 1:s
    subplot(t, s/t, i);
    imshow(SMatrix1(:,:,i), []);
    %title(sprintf('Attenuated Image %d', i));
    title(sprintf('Mosaic %d', i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These lines allow plotting each mosaic separately for visualization
% as independent images respecting their numbering and selection for
% focused study.
%for i = 1:s
%    A = SMatrix1(:,:,i);
%    figure;
%    imshow(A);
%    title(sprintf('Mosaic Num %d', i));
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main processing loop: Apply topological variational snake algorithm
% to each mosaic independently
% This extension of DetTACE allows batch processing of all image tiles
for i = 1 : s

    % Load image and detect maxima using topological approach
    % TACE_Max_NObj identifies candidate object centers based on
    % brightness and topological features in Sobolev spaces
    [maxima, enhanced_img] = TACE_Max_NObj(SMatrix1(:,:,i), n);

    % Plot the spots or ROI zones of each object
    %figure;
    %imshow(enhanced_img)

    % Segment objects using variational snake method
    % TACE_SegmentNObj applies the topological isophote method to
    % create binary masks that converge to object boundaries
    [BinaryMask, IndividualMasks] = TACE_SegmentNObj(SMatrix1(:,:,i), n, maxima);

    % Display binary mask for the current mosaic
    % BinaryMask contains the union of all individual object masks
    % This provides an overview of detected objects in the mosaic
    figure
    imshow(BinaryMask); 
    title(sprintf('Binary Mask from mosaic %d', i));

    % Visualization of centroids in the current mosaic
    % Overlay detected maxima (object centers) on the original mosaic
    figure
    imshow(SMatrix1(:,:,i), []); hold on;
    plot(maxima(:,2), maxima(:,1), 'r+', 'MarkerSize',15);
    title(sprintf('Original Image with Maxima - Mosaic %d', i));

    % Display individual masks for each detected object
    % This section is commented out but can be enabled for detailed
    % inspection of individual object segmentation
    %figure;
    %for i = 1:20
        %    imshow(IndividualMasks{i});
        %    title(sprintf('Binary mask object %d',i));
    %end

    % Overlay mask on original image
    % This image plots circular contours to which the variational method converges
    % The optimal boundary (from binary mask) and an extended boundary
    % (containing some noise) create a ring-like appearance
    % The reddish region represents the area between both edges calculated
    % variationally in TACE_Max_NObj and TACE_SegmentNObj
    figure;
    imshow(SMatrix1(:,:,i), []);
    hold on;
    visboundaries(BinaryMask, 'Color', 'r', 'LineWidth', 1);
    title(sprintf('Segmented Objects Overlay - Mosaic %d', i));

end