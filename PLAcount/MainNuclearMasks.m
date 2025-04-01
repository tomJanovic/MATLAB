close all
clear all
clc

%% Loading files
data = loadNucleiFiles(); %load Nuclei files
n = length(data(1,:)); %number of files loaded

%% Saving folder
pathname = [data{1,1} 'Nuclear masks/']; %pathname
if ~exist(pathname,'dir') %if folder does not exist
    mkdir(pathname); %create a folder
end

%% Processing files
for i = 1:n
    img = im2double(data{3,i}); %change to double
    
    %Normalization [0-1]
    imgNorm = img-min(min(img));
    imgNorm = imgNorm./max(max(imgNorm));
    
    %Filtering and blurring image
    imgBlurred = medfilt2(imgNorm); %median filtering
    blurMaskAverage = ones(15,15)/225; %average blur mask
    imgBlurred = conv2(imgBlurred, blurMaskAverage, 'same'); %convolution with an average blur mask
    
    %Segmenting nuclei
    mask = imbinarize(imgBlurred, 'global'); %thresholding with a global Otsu threshold
    mask = imfill(mask,'holes'); %filling holes
    mask = imopen(mask, strel('disk', 5)); %opening image
    mask = bwareaopen(mask, length(mask(:,1))); %select only regions with more pixels than the width of the image

    %Saving mask
    imwrite(mask, [pathname data{2, i}], 'tiff');
end