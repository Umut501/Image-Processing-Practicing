%% noisy1 recovering
% Read Image and check the image with imtool
    f = im2gray(imread('noisy1.png'));
    %imtool(f);
    
% Take a constant gray level block from the image and Plot the block 
    block = f(553:731,348:372);
    
% Get the histogram of the block 
    h = imhist(block);

% Design the gaussian pdf for uniform noise
    m = mean(double(block(:))); % mean
    v = std(double(block(:)));  % variance
    x = [0:255];

    ng = h;
%deciding values from histogram
    ng(ng<232)=0;
    ng(174:232)=71;
    ng(ng>169)=0;

% Figure the image histogram and uniform noise line    
    figure;
    subplot(1,3,1);
    bar(x,h);
    hold on;
    plot(x,ng,'r-','linewidth',2);
    title("histogram and uniform noise line");

% Image has additive noise type, that is Uniform Noise acording to the histogram
% Midpoint Filter - Works best for Gaussian and Uniform noise.
% I also implemented Arithmetic Mean filter to make image smoother

% padding image with zeros to apply the filter
paddedImage = [zeros(1,size(f,2)+2); zeros(size(f,1),1) f zeros(size(f,1),1) ;zeros(1,size(f,2)+2)];

% applying arithmetic mean filter to make image smoother
for j=1:size(paddedImage,1)-2%till reaching the end of the y dimension j is incremented
    for i=1:size(paddedImage,2)-2%till reaching the end of the x dimension i is incremented
        v(j,i) = sum(paddedImage(j:j+2,i:i+2),"all")./9;
     end
end

% applying midpoint filter to remove uniform noise
for j= 1:size( v,1)-2 %till reaching the end of the y dimension j is incremented
    for i=1:size( v,2)-2 %till reaching the end of the x dimension i is incremented
        maxP = max(max(v(j:j+2,i:i+2)));
        minP = min(min(v(j:j+2,i:i+2)));
        ImageOutput(j,i) = ((maxP+minP)./2);
     end
end

subplot(1,3,2);
imshow(uint8(ImageOutput),[]);
title("recovered image 1");
%Saving output image as a png file
imwrite(uint8(ImageOutput),"recovered1.png");

%Differences (subtraction) between edges of noisy and reconstructed images (subtraction).
paddedImageOutput = [zeros(1,size(ImageOutput,2)+2); zeros(size(ImageOutput,1),1) ImageOutput zeros(size(ImageOutput,1),1) ;zeros(1,size(ImageOutput,2)+2)];
dif = double(f) - double(paddedImageOutput);% subtracting to see edges

subplot(1,3,3);
imshow(dif,[]);%showing filtered image
title("Differences between edges image 1");

%% noisy2 recovering
%Read Image and check the image with imtool
    f = im2gray(imread('noisy2.png'));
    %imtool(f);

% We need to check the frequency domain for periodic noise
fimg = fft2(f);
cfimg = fftshift((fimg));% center image in frequency domain
%imtool(log(1+abs(cfimg)),[]);%check location of noise

%There are periodic noises in frequency domain we need to apply Band
%Reject Filter to remove them
%I applied Gaussian Band Reject Filter in the following function
M=512; N=512; Do=65; n=2;
[Hbp,Hbr] = BandPassAndRejectFilters('gaussian',M,N,Do,n,60);

%applying filter to image's frequency domain
cfimg = cfimg.*Hbr;
cfimg = double(cfimg);

figure;
subplot(1,3,1);
imshow(log(1+abs(cfimg)),[]);%show filtered frequency domain
title("filtered frequency domain");

u= real(ifft2(ifftshift(cfimg)));%convert image back to spatial

subplot(1,3,2);
imshow(uint8(u));%showing filtered image
title("recovered image 2");
%Saving output image as a png file
imwrite(uint8(u),"recovered2.png");

%Differences (subtraction) between edges of noisy and reconstructed images (subtraction).
dif = double(f)- double(u);
subplot(1,3,3);
imshow(dif,[]);%showing filtered image
title("Differences between edges image 2");


%% noisy3 recovering
%Read Image and check the image with imtool
    f = im2gray(imread('noisy3.tif'));
    %imtool(f);

% We need to check the frequency domain for sinusoidal noise
fimg = fft2(f);
cfimg = fftshift((fimg));% center in frequency domain
%imtool(log(1+abs(cfimg)),[]);%I checked the location of noise

%There is sinusodial noises in frequency domain we need to apply Notch Reject Filter
%I applied ideal Notch reject filter in the following function
M=246; N=168; Do=8; n=2;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,42,-28,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,40,28,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,-40,-28,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,-40,28,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,-80,-20,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,80,20,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,-80,20,n);
cfimg = cfimg.*Hnr;
[Hnp,Hnr] = NotchPassAndRejectFilters('ideal',M,N,Do,80,-20,n);
cfimg = cfimg.*Hnr;

cfimg = double(cfimg);
figure;
subplot(1,3,1);
imshow(log(1+abs(cfimg)),[]);%show filtered frequency domain
title("filtered frequency domain");

u= real(ifft2(ifftshift(cfimg)));%convert image back to spatial
subplot(1,3,2);
imshow(uint8(u));
title("recovered image 3");
%Saving output image as a png file
imwrite(uint8(u),"recovered3.png");

%Differences (subtraction) between edges of noisy and reconstructed images (subtraction).
dif = double(f)- double(u);
subplot(1,3,3);
imshow(dif,[]);%showing filtered image
title("Differences between edges image 3");



%% This function computes frequency domain bandpass and bandreject filters for gaussian Band Pass And Reject Filters
%   H = BandPassAndRejectFilters(TYPE, M, N, D0, n) creates a bandpass, Hbp, 
%   and a bandreject filter, Hbr, of the specified TYPE and size (M-by-N). 
function [Hbp,Hbr] = BandPassAndRejectFilters(TYPE,M,N,D0,n,W)
%   Valid values for TYPE, D0, and n are:
%   'gaussian'      Gaussian filter with cutoff frequency D0 and width W.

% Compute the distances D(U, V)
    for u=1:M
        for v=1:N
         D(u,v)=((u-(M/2))^2 + (v-(N/2))^2 )^(1/2);
         end
    end
    
% Compute the filter   
    switch TYPE
        case 'gaussian'
            Hbr=1-(exp(-((D.^2-D0.^2)./(D.*W)).^2));
            Hbp= 1-Hbr;
    end
end
%% This function computes frequency domain  filters for notch pass and reject filters 
function [Hnp,Hnr] = NotchPassAndRejectFilters(TYPE,M,N,D0,uk,vk,n)
%Computes frequency domain notchpass and notchreject filters.
%   H = NotchPassAndRejectFilters(TYPE, M, N, D0, n) creates a notchpass,  
%   Hnp,and a notchreject filter, Hnr, of the specified TYPE and 
%   size (M-by-N). 
%
%   Valid values for TYPE, D0, and n are:
%   'ideal'         Ideal  filter with cutoff frequency D0 and centers(uk,vk).
%
%   'butterworth'   Butterworth filter of order n, cutoff D0 and centers(uk,vk).
%
%   'gaussian'      Gaussian filter with cutoff frequency DO and centers(uk,vk).

% Define D(u,v)
  
    for u=1:M
        for v=1:N
         Dkp(u,v)=((u-(M/2)-uk)^2 + (v-(N/2)-vk)^2 )^(1/2);
         Dkn(u,v)=((u-(M/2)+uk)^2 + (v-(N/2)+vk)^2 )^(1/2);
         end
    end
    
% Compute the filter   
    switch TYPE
        case 'ideal' 
            Hnp=zeros(M,N);
            Hnp(Dkp<=D0)=1;
            Hnr=1-Hnp;
        case 'gaussian'
         Hnr = 1-exp((-1/2).*((Dkn.*Dkp)./(D0.*D0)));
         Hnp = 1-Hnr;
    end
end
