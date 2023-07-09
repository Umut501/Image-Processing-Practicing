%I converted image to a double to make operations on it 
%Saved it into Image2 variable
Image2 = double(imread('Image2.png'));

%Checking converted Image with imtool to decide noise type
imtool(Image2);
%Image has salt and pepper noise to remove it needed filter type is median filter

% padding Image2 with zeros to apply the filter
paddedImage2 = [zeros(1,size(Image2,2)+2); zeros(size(Image2,1),1) Image2 zeros(size(Image2,1),1) ;zeros(1,size(Image2,2)+2)];

% applying median filter
for j= 1:size(Image2,1)%till reaching the end of the y dimension j is incremented
    for i=1:size(Image2,2)%till reaching the end of the x dimension i is incremented
        %taking the median of the 3 by 3 piece and putting it to output image
        %while doing this each pixel value is converted to uint8 since they
        %are converted to double at the beginning of the process
        Image2Output(j,i) = uint8((median(paddedImage2(j:j+2,i:i+2),"all"))); 
        Image2Output_(j,i) = (median(paddedImage2(j:j+2,i:i+2),"all"));%this output is created to use in the next part of the question it is not converted to uint8 for the purpose of usage  
     end
end

%showing the output image
figure;
subplot(2,2,1);
imshow(imread('Image2.png'));
title("Image2");
subplot(2,2,2);
imshow(Image2Output);
title("Image2Output");

%saving output image as a png file
imwrite(Image2Output,"Image2Output.png");

%-----------------------------------------------------------------------------------------

%For the edge detection the needed filter types are Horizontal Sobel and Vertical Sobel filters

%Creating Horizontal Sobel and Vertical Sobel filters
Horizontal_Sobel_Filter = [-1 -2 -1; 0 0 0; 1 2 1];
Vertical_Sobel_Filter = [-1 0 1;-2 0 2; -1 0 1];

% Padding Image2Output_ with zeros to apply the filter
% Image2Output_ created in the first part of the question
paddedImage2Output = [zeros(1,size(Image2,2)+2); zeros(size(Image2,1),1) Image2Output_ zeros(size(Image2,1),1) ;zeros(1,size(Image2,2)+2)];

% applying sobel filters
for j= 1:size(Image2,1)%till reaching the end of the y dimension j is incremented
    for i=1:size(Image2,2)%till reaching the end of the x dimension i is incremented

        %To apply filters I multiplied them with the given outputs and take the sum of created array. 
        %It is equivalent operation of the correlation
        %After that to create final image I added their absolute values as the formula states in the lecture notes

        %Both the original image and the image without salt and pepper noise is processed with sobel filters 
        %And their outputs are recorded in uint8 format

        Image2Out_sobel_h(j,i) = sum(paddedImage2Output(j:j+2,i:i+2).*Horizontal_Sobel_Filter(1:3,1:3),"all");
        Image2Out_sobel_v(j,i) = sum(paddedImage2Output(j:j+2,i:i+2).*Vertical_Sobel_Filter(1:3,1:3),"all");
       
        Image2Output_Edges(j,i) = uint8(abs(Image2Out_sobel_h(j,i))+abs(Image2Out_sobel_v(j,i)));
        
        Image2_sobel_h(j,i) = sum(paddedImage2(j:j+2,i:i+2).*Horizontal_Sobel_Filter(1:3,1:3),"all");
        Image2_sobel_v(j,i) = sum(paddedImage2(j:j+2,i:i+2).*Vertical_Sobel_Filter(1:3,1:3),"all");

        Image2_Edges(j,i) = uint8(abs(Image2_sobel_h(j,i))+abs(Image2_sobel_v(j,i)));

    end
end

%the edges of the image without salt and pepper noise is shown here
subplot(2,2,4);
imshow(Image2Output_Edges);
title("Edges of Image2Output");

%the edges of the original image is shown here
subplot(2,2,3);
imshow(Image2_Edges);
title("Edges of Image2");
