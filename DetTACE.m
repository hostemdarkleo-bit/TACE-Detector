%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled "Semi-Automatic model based on topologic 
%    isophotes in Sobolev spaces to identify objects in astrophysical images:
%    Base of an AI algorithm". 
%
%    One of the main challenges of our work is to compete in some
%    routines with SExtractor et al. Bertin (10.1051/aas:1996164) for
%    object identification and segmentation by optimizing contour
%    construction. Based on closed curves in Sobolev spaces to apply
%    variational analysis and build binary masks that converge to the
%    outer boundary of objects and thus segment them by identifying
%    their centroid.
%
%    Our proposal is modular in both Matlab and Python by dividing it
%    into functions that must be loaded before executing the main routine
%    called "DetTACE.m". Before executing this routine, the following
%    functions must be loaded and executed (to ensure they are active):
%    IntFilImg.m, LogFilterIm.m, SMosaics.m, create_mask.m,
%    estimate_radius.m, TACE_Max_NObj.m, TACE_SegmentNObj.m.
%
%    Corresponding Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%    Dra. Bolivia Teresa Cuevas Otahola
%    bcuevas@astro.unam.mx;          b.cuevas.otahola@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read metadata
Img = fitsread('C:\MathematicsNew\Bolivia.fits');
%Img = fitsread('C:\MathematicsNew\JWebDataPaper\Test Images\h_m82_i_s05_drz_sci.fits');
%Img = img.data;

figure;
imshow(Img)
title('Original study Image');

% Attenuated image considering only its intensity, with attenuation factor "AtFac = 12"
AtFac = 5;
Img = Img/AtFac;

% Partition into 6 tiles (s=6), arranged in 2 rows (t=2), meaning we create
% "s" divisions (converting the original image into a mosaic of s
% sub-images) for "t" number of sub-image rows in the mosaic
s = 9;
t = 3;

% Number of objects of interest based on brightness per sub-image in the mosaic
n = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partition where mosaics are aligned and arranged in descending order
% from the first entry of the top row to the last row from beginning to end
[SMatrix1, MtPxData1] = SMosaics(Img, s, t);

% Visualize tile positions
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

% Test filters on a central sub-image with higher saturation
% Dynamic filters are tested where in the bottom bar we can appreciate
% the filter dynamics and the filter result. When closing the window,
% this result is saved in the variable

% Selected mosaic to be analyzed and standardized by applying
% "asinhStretchInteractive" analysis to the entire image, as a deeper
% processing with optimal results.
image = SMatrix1(:,:,5);       % Mosaic with highest saturation
% If we replace the previous mosaic with the following image segment, we 
% can obtain the results in Fig. 1 (b)-(d)
%image = imcrop(Img, [8544.5 5522.5 643 586]);

Image = im2double(image);      % Convert mosaic to double [0,1]

% Simple intensity filter
% We can choose a simple intensity filter based on Fourier Transform
% and subsequently a logarithmic filter. Choosing either one filter
% or both on the test image.
AtImage = IntFilImg(Image);
AtImage0 = LogFilterIm(AtImage);

% If we wish to apply these filters later to the analysis mosaic
% to discriminate more objects in the following steps, we apply them
% with the respective functions "IntFilImg" and "LogFilterIm".
% Then we choose by our conviction and needs which image to work with
% in the variable "image".
image = AtImage;

figure;
imshow(image);
title('Test Image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load image and detect maxima
%[maxima, ~] = TACE_Max_NObj(image, n);
[maxima, enhanced_img] = TACE_Max_NObj(image, n);

% Plot the spots or ROI zones of each object
figure;
imshow(enhanced_img)

% Segment objects
[BinaryMask, IndividualMasks] = TACE_SegmentNObj(image, n, maxima);

% Plot masks of each object in the entire mosaic. We obtain each mask
% per mosaic, but that depends on the user if they want or need to
% focus on a single object per ROI
imshow(BinaryMask)

% Visualization of centroids per mosaic
figure
imshow(image, []); hold on;
plot(maxima(:,2), maxima(:,1), 'r+', 'MarkerSize',15);
title('Original Image with Maxima');

% Visualization of centroids per mosaic and of the mask, in both
% cases for the total objects
figure;
subplot(1,2,1);
imshow(image, []); hold on;
plot(maxima(:,2), maxima(:,1), 'r+', 'MarkerSize',15);
title('Original Image with Maxima');
subplot(1,2,2);
imshow(BinaryMask);
title(sprintf('Segmented %d Objects',n));

% Show individual masks
figure;
for i = 1:min(n,4)
    subplot(2,2,i);
    imshow(IndividualMasks{i});
    title(sprintf('Object %d',i));
    hold on;
    plot(maxima(i,2), maxima(i,1), 'r+', 'MarkerSize',10);
end

% Overlay mask on original image. This image plots circular contours
% to which we converge, considering the optimum taken in the binary mask
% and a more extensive or larger one containing some noise.
% Therefore, it appears as a ring since the reddish region is an area
% drawn between both edges calculated variationally in TACE_Max_NObj
% and TACE_SegmentNObj
figure;
imshow(image, []);
hold on;
visboundaries(BinaryMask, 'Color', 'r', 'LineWidth', 1);
title(sprintf('Segmented %d Objects Overlay', n));