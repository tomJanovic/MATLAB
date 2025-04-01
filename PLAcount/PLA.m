close all
clear all
clc

%% Selecting the DAPI and PLA channels
% 1 - RED
% 2 - GREEN
% 3 - BLUE

PLA = 1; %PLA signal is by default RED (1) channel
DAPI = 3; %DAPI signal is by default BLUE (3) channel

%% Loading files
data = loadPLAfiles(); %load files
PLAdata = struct([]); %allocate struct for output data
vecPLA = []; %allocate the vector of PLA foci

if iscell(data) == 0 %if only one file is loaded
    n = 1;
else
    n = length(data); %if more files are loaded
end

blurMaskGauss = [1 2 4 2 1
                 2 4 8 4 2
                 4 8 16 8 4
                 2 4 8 4 2
                 1 2 4 2 1]/100; %gaussian blur mask
             
blurMaskAverage = ones(5,5)/25; %average blur mask
               
for i = 1:n
    if n == 1
        imgPLA = im2double(data); %get image
    else
        imgPLA = im2double(data{i}); %get image
    end
    REDch = imgPLA(:,:,PLA); %pick PLA channel
    BLUEch = imgPLA(:,:,DAPI); %pick DAPI channel
    
%% Adding artificial blur
%     REDch = imnoise(REDch,'gaussian', 0.01);
%     BLUEch = imnoise(BLUEch,'gaussian', 0.01);
    
%% Gaussian blur
    REDch = medfilt2(REDch);
    REDch = conv2(REDch, blurMaskGauss, 'same'); %convolution with gaussian blur mask
    BLUEch = conv2(BLUEch, blurMaskAverage, 'same'); %convolution with average blur mask
    
%% Normalization [0-1]
    REDch = REDch-min(min(REDch));
    REDch = REDch./max(max(REDch));
    BLUEch = BLUEch-min(min(BLUEch));
    BLUEch = BLUEch./max(max(BLUEch));

%% Obtaining DAPI mask
    otsuDapiThresh = graythresh(BLUEch)/2; %global threshold as a half of optimal Otsu threshold
    maskDAPI = im2bw(BLUEch, otsuDapiThresh); %thresholding
    maskDAPI = imfill(maskDAPI,'holes'); %filling holes
    maskDAPI = imopen(maskDAPI, strel('disk', 5)); %opening image
    maskDAPI = imclearborder(maskDAPI); %removing nuclei that intersect with the image edge
    maskDAPI = bwareaopen(maskDAPI, length(maskDAPI(:,1))); %select only regions with more pixels than width of the image

%% Obtaining PLA signal only from nucleus
    PLAsignal = REDch.*maskDAPI; %multiply the DAPI mask with RED channel
    plaThresh = mean(mean(PLAsignal(maskDAPI))) + graythresh(PLAsignal(maskDAPI)); %mean value + otsu of PLA signal within the nuclei
    PLAsignal = REDch.*(im2bw(PLAsignal, plaThresh)); %subtract the background with plaThresh
%     PLAsignal = imopen(PLAsignal, strel('disk', 1)); %opening image
    localMaxima = imregionalmax(PLAsignal, 8); %obtain local maxima

%% Watershed segmentation
    inverseGradientMap = imhmin(1-bwdist(~maskDAPI), 6); %distance transform of binary image (negative nucleus) and supressing of regional minima
    segmentedRegions = watershed(inverseGradientMap); %watershed segmentation

%% Count PLA signals for each nucleus
    nucleiCount = max(max(segmentedRegions)); %count number of regions in segmented image
    totalCountPLA = zeros(nucleiCount,1); %allocate PLA count
    
    for j = 1:nucleiCount %for each nucleus
        selectedRegion = zeros(size(maskDAPI)); %allocate region
        selectedRegion(segmentedRegions == j) = 1; %pick only one region
        selectedRegion = (selectedRegion + maskDAPI) - 1; %add maskDAPI and subtract 1
        selectedRegion(selectedRegion < 0) = 0; %each pixel below 0 assign to 0
        PLAfoci = selectedRegion.*localMaxima; %multiply the region with localMaxima to obtain only PLA foci from one nucleus
        PLACount = bwconncomp(PLAfoci);
        totalCountPLA(j,1) = PLACount.NumObjects; %store the number
    end
    PLAdata{i} = totalCountPLA; %store the amount of PLA foci to the database
    vecPLA = [vecPLA totalCountPLA'];
end

%% Statistics
AveragePLAforNuclei = double(sum(vecPLA))/double(length(vecPLA)); %count average number of foci per nucleus

%% Overlay of DAPI and PLA signal with segmented regions
nucleusSegmentedImage = zeros(size(imgPLA));
nucleusSegmentedImage(:,:,3) = maskDAPI.*BLUEch;
nucleusSegmentedImage(:,:,2) = segmentedRegions == 0;
nucleusSegmentedImage(:,:,2) = nucleusSegmentedImage(:,:,2) + localMaxima;
nucleusSegmentedImage(:,:,1) = REDch;

%% Figures
figure(1) %show segmented nuclei
title('Segmented nuclei regions');
imshow(nucleusSegmentedImage, [])

figure(2) %show original PLA signal and detected PLA foci
subplot(1,2,1)
imshow(imgPLA(:,:,1), [])
title('Original PLA signal');
subplot(1,2,2)
imshow(localMaxima)
title('Detected local maxima of PLA signal');

figure(3) %show boxplot analysis
boxplot(vecPLA)
title('PLA in U251 cell line')
ylabel('Number of PLA foci')

%% Display statistics
disp(['Amount of cells: ',num2str(length(vecPLA))])
disp(['Average of PLA foci per nucleus: ',num2str(AveragePLAforNuclei)])
disp(['Median of PLA foci per nucleus: ',num2str(median(vecPLA))])

