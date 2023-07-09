%Saving Image1.png into the image1 variable
%After converting it to a grayscale image
%I have done this to ensure the pixel values are between 0 and 255.
image1 =im2gray(imread("Image1.png"));

%showing the Image1 with its historam using figure,imshow and subplot commands
figure;
subplot(1,2,1);% putting histogram of the image at the first place of figure
imhist(image1);%creates a histogram for the image1
title("Image1 Histogram");
subplot(1,2,2);% putting the image1 at the second place of the figure
imshow(image1);
title("Image1");

%After checking the histogram data, I saw it has low-contrast values
%I decided to use histogram equation since it enhances the contrast of the images as mentioned in the lecture notes:
%"... increases the dynamic range of the gray-levels in a low-contrast image to cover full range of gray-levels"
Image1Output = histeq(image1);%applying histogram equation to Imge1

%Saving output image as a png file
imwrite(Image1Output,"Image1Output.png");

%showing the Image1Output with its historam using figure,imshow and subplot commands
image1Output = imread("Image1Output.png");
figure;
subplot(1,2,1);
imhist(Image1Output);%creates a histogram for the image1
title("Image1Output Histogram");
subplot(1,2,2);% putting the Image1Output at the second place of the figure
imshow(Image1Output);
title("Image1Output");
