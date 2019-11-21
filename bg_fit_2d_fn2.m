function [bg_surf]=bg_fit_2d_fn2(img_of_interest,x0,y0,r0)

%Description: function to generate a background surface based on what hits 
%the camera outside the NA. This surface is then subtracted from the BFP 
%image. This method is an improvement over previous methods, including 
%subtracting a bg 1) val based on a patch of the image in a single location
%2) val based on multiple patch locations or 3) clumsy attemts to make a 
%surface by linaer interpolation. 

%INPUTS: the image, as well as origin (x,y) and radius of the circle to
%exclude from the background surface determinat, usually a cirlce a bit
%bigger than the NA to squash bright stuff at the edge

%OUTPUTS: background surface to subtract

%zjs, cleaned up 10-21-2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=img_of_interest;
mask=ones(size(M));
[h,w] = size(M);
[X,Y] = meshgrid(1:w,1:h); %creates a grid where the vals are pix locations

%debug image in plus overlay of circle, surf derived from what's not inside
%{
figure(17)
imagesc(M)
caxis([0 25])
hold on
circle(x0,y0,r0)
hold off
%}

%loop over pixels, exclude if they fall in circle. 
for n=1:h
    for m=1:w        
    r_squared=(m-x0)^2+(n-y0)^2;
    if r_squared<r0^2
       mask(n,m)=nan; 
       M(n,m)=nan;%make inside the circle nans
       X(n,m)=nan;
       Y(n,m)=nan;
    end
    
    end
end

%debug images, masked input, masked grids
%{
figure(12)
imagesc(M)
drawnow
figure(13)
imagesc(X)
figure(14)
imagesc(Y)
%}

%hmm, organize images and grids to vectors
X=X(:);
Y=Y(:);
Mv=M(:);
%discard NANs
X=X(~isnan(X));
Y=Y(~isnan(Y));
Mv=Mv(~isnan(Mv));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting happens here, fit to vectorized data, note: there are different 
%options for fitting functions, settled on something relatively flat, poly11
%sf = fit([X, Y],Mv,'poly23')
%sf = fit([X, Y],Mv,'poly33')
%sf = fit([X, Y],Mv,'poly22')
sf = fit([X, Y],Mv,'poly11');
    
%debug fig- tries to plot surf and bg points at the same time, doesn't work
%very well:
%figure(15)
%plot(sf,[X,Y],Mv)
%drawnow

%reconstitue a grid for the whole image and generate the surface from it: 
[X,Y] = meshgrid(1:w,1:h);
M_fit=sf(X,Y);

%debug image, image of the surface: 
%figure(16)
%imagesc(M_fit)
%M=M';

bg_surf=M_fit;
end