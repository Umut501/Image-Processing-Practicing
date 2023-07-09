%Saving Image3.png into the Image variable
Image = imread('Image3.png');
%I decided to use treshold to spot the birghthest star since it allows us
%to segment image and understand it better
%Checking Image with imtool to decide a treshold value for segmentation
%I decided treshold value as 126
imtool(Image);

%I converted Image to a double variable to make operations on it
Image3 = double(Image);

%Converting pixels to pure black if they are less than the treshold value
Image3(Image3<126)=0;

%Padding Image3 with zeros to apply a filter
paddedImage3 = [zeros(1,size(Image3,2)+2); zeros(size(Image3,1),1) Image3 zeros(size(Image3,1),1) ;zeros(1,size(Image3,2)+2)];

% Applying median filter to remove little stars, since they look like salt noise
for j= 1:size(Image3,1)%till reaching the end of the y dimension j is incremented
    for i=1:size(Image3,2)%till reaching the end of the x dimension i is incremented
        %taking the median of the 3 by 3 piece and putting it to output image
        %while doing this each pixel value is converted to uint8 since they
        %are converted to double at the beginning of the process
        Image3Output(j,i) = uint8(median(paddedImage3(j:j+2,i:i+2),"all"));
    end
end

%showing the output image
figure;
subplot(1,2,1);
imshow(Image)
title("Image3");

subplot(1,2,2);
imshow(Image3Output);
title("Image3Output");

%saving output image as a png file
imwrite(Image3Output,"Image3Output.png");
