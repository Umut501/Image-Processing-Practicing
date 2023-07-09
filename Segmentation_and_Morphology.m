function A3_2385201(image)
    % Reading the original image
    original_image = imread(image);

    % Segmenting the image
    segmented_image = segmentation(original_image);

    % Counting the number of big and small half-cut oranges
    counting(segmented_image);
end

%% Median filtering and Otsu's method are used in this function to segment the oranges here because 
% median filtering smooths the image and removes the noise, 
% and Otsu's method finds the optimal threshold value to separate the foreground and background pixels.

function segmented_image = segmentation(original_image)
    % Reading the image
    f = rgb2gray(original_image);

    % Pre-processing the image using median filtering
    f = medfilt2(f);

    % Using Otsu's method (graythresh built-in function) to get a threshold 
    % value between 0 and 1
    thresh = graythresh(f);

    % Converting the otsu threshold value to be between 0-255
    thresh = thresh.*255;

    % Applying bi-level thresholding (See slide 35) and create new image sb
    sb=f;
    sb(sb<thresh)= 0;
    sb(sb>=thresh)= 255;

    % Creating Structuring Element S with 'disk' and size 3
    S = strel('disk', 3);

    % Applying Opening = Erosion and Dilation
        % Erode image f using imerode built-in function and Structuring 
        % Element S and obtain image E1
        E1 =imerode(sb, S);
        % Dilate image E1 using imdilate built-in function and Structuring 
        % Element S and obtain image D1
        D1 = imdilate(E1, S);

     % Inverting the image E2 since bwconncomp checks white areas
        E2 =255 - D1;

     % Fill in any incomplete regions in the image using imfill
        E2 = imfill(E2, "holes");

    % Showing the all of the steps together
    figure;
    subplot(2, 2, 1);
    imshow(original_image); title('Original Image');
    subplot(2, 2, 2);
    imshow(f); title('Grayscale and Median Filtered Image');
    subplot(2, 2, 3);
    imshow(255-sb); title('Binary Image');
    subplot(2, 2, 4);
    imshow(E2); title('Applying Opening and Filling');

    % Capturing the image data from the figure
    frame = getframe(1); % gets the firstly displayed frame

    % Saving the image as a PNG using imwrite
    imwrite(frame.cdata, 'Figure1.png', 'png'); 

    % Returning the segmented image
    segmented_image = E2;
end

%% I use bwconncomp function here which finds and counts the connected components in a binary image
% and cellfun function to find how many pixels exist in a component.
% After that I used it to find a number of small oranges and big oranges
% separately by calculating the avrage number of pixels in a component.

function [num_big, num_small] = counting(segmented_image)
    % Finding the number of connected components in the segmented_image using
    % bwconncomp built-in function
    C =  bwconncomp(medfilt2(segmented_image)); %I also apply median filter here to remove if there is any noise caused by the small components

    % Geting the value of connected components
    numComponents = length(C.PixelIdxList);

    % Geting the number of pixels in each connected component  
    numPixels = cellfun(@numel,C.PixelIdxList);

    % Geting the largest number of pixels of the connected components   
    [biggest,~] = max(numPixels);

    % Geting the minimum number of pixels of the connected components   
    [smallest,~] = min(numPixels);

    % Geting the average number of pixels of the connected components   
    average = (smallest+biggest)/2;

% Creating counters for oranges
big_orange_counter = 0;
small_orange_counter = 0;

% Checking all connected components via this loop
for i = 1:numComponents

% Getting the number of pixels in the current connected component
numPixels = size(C.PixelIdxList{i}, 1);

 % Comparing the size of the avarage connected component with the biggest connected component
  if numPixels > average
      % Incrementing the counter for big oranges
      big_orange_counter = big_orange_counter + 1;
  end

  if numPixels <= average
      % Incrementing the counter for small oranges
      small_orange_counter = small_orange_counter + 1;
  end
end

% Creating a figure showing the segmented image and the number of big and small oranges
figure;
imshow(segmented_image);
% Here I use num2str function to convert numbers to string to display them in the title
title(['Number of big oranges: ' num2str(big_orange_counter) ', number of small oranges: ' num2str(small_orange_counter)]);
   
% Capturing the image data from the figure
frame = getframe(2); % gets the secondly displayed frame

% Saving the image data to a PNG using imwrite
imwrite(frame.cdata, 'Figure2.png', 'png');
end