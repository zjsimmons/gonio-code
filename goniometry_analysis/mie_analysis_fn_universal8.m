function [lists_out,mie_center,mask_out,theta_vector,radial_avg_downsampled_matrix] =...
    mie_analysis_fn_universal8(img,bg,mask_in,mie_x0,mie_y0,ancillary_inputs,analysis_settings)

%disp('a')
%tic

%DESCRIPTION: This is the main analysis program for analyzing BFP images. 

%This program is pretty complicated, the general idea is the following:
%The program requires a BFP image, background image as well as other input 
%parameters: the pupil plane to angle conversion, optional mask of areas to
%discard from analysis, the center of the forward-scattered pattern as well
%as other settings such as the location and size of the BFP.

%The BFP image is mapped to angle space via asin, i.e. 'undoing' the
%distortion of the objective. The angle space is rotated via Euler
%rotations so that all transformed images have the same orientation. The
%program outputs a variety of diagnostic outputs, as well as center
%locations, and mask that can be fed back in to expedite re-running the
%program. 

%Originally, the main output was a trace of the average over phi. Later
%functionality to output a binned version with both theta and phi was
%incorporated. 

%Subsequent versions also included logic/effort was invested in optimizing
%the program to get it to run as fast as possible, to facilitate analyzing
%lots of images faster. 

%INPUTS: 
    %img - input image to analyze
    %bg - background image, can also be zero
    %pp_coord,angle_out normalized pupil coordinates and corresponding
    %angles in sample (in H2O)
    %mask_in - for later runs, this lets you use alread generated 
    %user-chosen masks 
    %mie_x0,mie_y0 - pattern center, once found can be fed back in for
    %differerent runs, e.g. designate w/ set where centers are very clear
    %analysis_settings - number of user defined settings including: 

    %analysis_settings.analysis_circle_x - where to analyze circle
    %analysis_settings.analysis_circle_y
    %analysis_settings.analysis_circle_NA
    
    %analysis_settings.mask_circle_x - where to mask (can differ)
    %analysis_settings.mask_circle_y
    %analysis_settings.mask_r
    
    %analysis_settings.fw_on=1;%1=fw, 0=bw forward or backward analysis
    %analysis_settings.color=525; %for accompanying mie trace
    %analysis_settings.use_input_mask_on=1; %use mask input or user select
    %analysis_settings.center_lookup_on=1; %use input centers or user
    %analysis_settings.how_many_phi_bins=1; %for g(phi), fragile
    %analysis_settings.figures_off=1; %to turn figures off when doing a lot
    %of processing
    
%user inputs: 1) calibration inputs: NA circle, mask circle locations and
%sizes. 2) user-selected scattering parameter centers. 3)
%user-selected/generated masks for what of the BFP image to include in
%analysis. -all these user inputs are hammered out on the initial runs,
%when satisfied, they're fed back in. 

%OUTPUTS: 
    %theta_vector - trace averaged over phi- theta
    %radial_avg_downsampled_matrix - trace averaged over phi
    
    %mie_center_x - user selected pattern center and mask
    %mie_center_y
    %mask_out
    
    %x_todisp_ds - for 3D display, fragile 
    %y_todisp_ds
    %z_todisp_ds
    %r_todisp_ds
    
    %lists_out - phase function portion as fn of theta AND phi
    
    %g_avg - g averaged around phi

%FUNCTIONS CALLED: 
%[fList,pList] = matlab.codetools.requiredFilesAndProducts('mie_analysis_fn_universal5');
%fList'
%mie_analysis_fn_universal5.m:

%bg_fit_2d_fn2.m - generates the bg surf to subtract off the image 
%bin_it3c.m -takes the pixel brightness list and theta list and bins them
%circle.m - draws a circle
%circle_with_color.m - draws a circle with color option
%g_integral.m - calculates g frome the forward partial trace
%lininterp1.m - interpolates when mapping pixel to angle

%REVISION HISTORY (broad):
%mie_analysis_fn callable function to look at images
% originally written fall 2015 zjs

%1) this version has a plot option 
%2) also has inputs to define the circular field of view to include.. in
%analysis
%3) update end of dec 2015: incorporates off-center correction
%4) jan 2015: improved off-center correction
%5) this version is set to look in the backwards direction, NA, etc.

%6) need functionality for normalized bpp rather than pixel look-up

%3-30-2016, just tried cos(\theta) angle brightness compression correction,
%does the same thing as the functionality I had, that's comforting. Fully
%fixed 6-10-2016 

%4-28-2016, lets add logic to discard saturated pixels..
%perhaps if a fraction of saturated pixels at a given angle is too large
%they get discarded..

%July 2016, big fix, moved from shifting to rotation to correct for
%off-center, seems to fix the problem. 

%7-19-2016, added g(phi) functionality..

%Aug 2016, implemented better background subtraction functionality. Now
%generate a surface based on what is outside the NA and subtract that
%rather than 1) a value from a patch in a corner, 2) avg of 4 corner
%patches, 3)crude surface generated by averaging linear interpolations
%between patches. New surface method seems more powerful, cleaner as it's
%in a fn. 

%Oct 2016, tightened up the code, removed a lot of legacy and other garbage

%spring 2017, revised to add *real* 4pi output rather than averaging over
%phi

%zjs, cleaned up 8-15-2017

%big streamline, 2018.1.20, zjs
%-Improved binning functionality so it runs a lot faster. 
%-Removed a lot of legacy code

%Notes/To Do: 

%-Fix angle brightness compression to account for different media n?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%user defined: 
%many parameters are contained in analysis_settings structure array:
cir_x=analysis_settings.analysis_circle_x;
cir_y=analysis_settings.analysis_circle_y;
cir_rad_NA=analysis_settings.analysis_circle_NA;
mask_x=analysis_settings.mask_circle_x;
mask_y=analysis_settings.mask_circle_y;
mask_r=analysis_settings.mask_r;
fw_on=analysis_settings.fw_on;

hard_coded_on=analysis_settings.use_input_mask_on;%this turns on and off 
%selection mask and center finding. Note, center lookup must also be on
center_lookup_on=analysis_settings.center_lookup_on;
color=analysis_settings.color;
figs_off=analysis_settings.figures_off;
speed_mode=analysis_settings.speed_mode;

%how_many_phi_bins=1;
how_many_phi_bins=analysis_settings.how_many_phi_bins;

%also in ancillary_inputs 
pp_coord=ancillary_inputs.pp_coord;
angle_out=ancillary_inputs.angle_out;
thetas=ancillary_inputs.thetas;
phis=ancillary_inputs.phis;

%hard-coded:
img=double(img);
bg=double(bg);
%use less than the circle selected: 
%cir_rad=cir_rad_NA-50;
cir_rad=cir_rad_NA;

% use a smaller circle (than the NA) to define the edge of the mask?
use_mask_circle=1;
if use_mask_circle==1
    cir_rad_to_mask=mask_r;
else
    cir_rad_to_mask=cir_rad_NA;
end

% this let's you subtract background image.. if the input bg image is zero,
% this doesn't do anything otherwise:
img_for_sat_mask=img;

use_bg=1;
if length(bg)>1 && use_bg==1
  %  disp('subtract off bg') %for debugging
  
 % class(img)
 % class(bg)
  
   img=img-bg; 
end

%create sat mask for use later to discard sat pixels in bin_it fn.
sat_mask=zeros(size(img_for_sat_mask));
%sat_val=1023; 
sat_val=max(max(img_for_sat_mask));
%sat_val=2^16-1;
sat_mask(img_for_sat_mask==sat_val)=1;
%sat in either? to drive down big subtractions..
sat_mask(bg==sat_val)=1;

%figure(68)
%imagesc(sat_mask)
%drawnow

%this turns on and off center-finding logic for fine-tuning fiding the Mie
%center pattern. (consider removing)
do_looping=0;

%the water-glass circle..
    water_NA=cir_rad_NA*1.33/1.4;
    air_NA=cir_rad_NA*1/1.4;
    water_edge_turn_on=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Image Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bg subtraction, circle mask generation, selection of analysis region

%generate the background surface to subtract, based on stuff outside NA. 
%we have an offset here so that it's not contaminated by bright stuff on 
%the edge of the NA [NOTE: do we want this functionality still?]
if analysis_settings.subtract_bg_surf==1
    bg_img_surf=bg_fit_2d_fn2(img,cir_x,cir_y,cir_rad_NA+50);
else
    bg_img_surf=0;
end

%if you want constant or zero, can hard-code here: 
%bg_img_surf=0;    

img=img-bg_img_surf; 
img(img<0)=0;%for negative vals to zero

%discard stuff outside our fov..
[no_rows,no_columns]=size(img);
zero_mask=ones(size(img));
nans_mask=ones(size(img));
%toc %a

%disp('b')
%tic
%create nan and zero masks based on input mask circle
for n=1:no_rows
    for m=1:no_columns
        %y_fov=n-cir_y;
        %x_fov=m-cir_x;
        y_fov=n-mask_y;%these should coincide, probably doesn't matter how this is done..
        x_fov=m-mask_x;
        %create a mask based on the FOV circle of pixels to include:
        argument2=x_fov^2+y_fov^2;
        %   if argument2>(cir_rad)^2
        if argument2>(cir_rad_to_mask)^2
            %binary 1 or 0 mask
            zero_mask(n,m)=0;
            nans_mask(n,m)=nan;%this creates a mask that kills everything outside the mask circle...
        end
    end
end
%toc %b
    %to discard%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %kill pixels outside of field of view..
    %img=img.*zero_mask;

    %scale to fix angle compression brightness:
    %old angle compression correction...
    %img=img./image_scaling_mask;

    %scale image between 0 and 1:, don't want to do this- screws up stuff when
    %you're outside the NA circle..
    %img=(img-min(min(img)))/(max(max(img))-min(min(img)));


    %compare the angle compression in image..
    %{
    img2=img./image_scaling_mask;
    img2=(img2-min(min(img2)))/(max(max(img2))-min(min(img2)));
    figure(11)
    imagesc(log(img))
    figure(12)
    imagesc(log(img2))
    %}

%plot our image and add some stuff to look at, selection logic..
    %{
figure(30)
imagesc(log(double(img)));
hold on
circle(mask_x,mask_y,mask_r)
hold off
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check out our BFP image...

%disp('c - only takes time if displaying')
%tic

if figs_off~=1
    figure(2)
    hold on
    set(gca,'YDir','reverse') %seem to have to do this to get proper display
    g=imagesc(log(double(img)));
    circle(cir_x,cir_y,cir_rad_NA)
    plot(cir_x,cir_y,'bX')
    drawnow
    if use_mask_circle==1
        figure(2)
        circle(mask_x,mask_y,mask_r)%additional mask circle
    end
    %water-glass edge..
    if water_edge_turn_on==1
        circle_with_color(cir_x,cir_y,water_NA,'-b')
        circle_with_color(cir_x,cir_y,air_NA,'-b')
    end
    %colormap(gray)
    axis tight %seem to have to do this to get it to display properly
    set(gca,'YDir','reverse') %seem to have to do this to get proper display
end

%User input Mask selection...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hard_coded_on~=1
    disp('Select where you want to analyze')
    h=imfreehand;
    %position = wait(h); %add to have to click
    %this is a mask generated from the initially selected region:
    mask=createMask(h,g);
    nans_mask2=ones(size(img));
    nans_mask2(mask==0)=nan;
    
    %this is the nans mask that is initial nans mask multiplied by the selected
    %nans mask:
    nans_mask2=nans_mask.*nans_mask2;
    %if you want to look at the masks:
    %figure(7)
    %imagesc(nans_mask2)
    
    %exclude some other spots for example?:
    exclude_flag='y';
    prompt = 'Select regions to exclude? y/n ';
    exclude_flag = input(prompt,'s');
    while exclude_flag =='y'
        nans_mask_swap=ones(size(img));
        figure(2)
        %create a new region and mask, now the selected region is excluded.
        h2=imfreehand;
        mask2=createMask(h2,g);
        nans_mask_swap(mask2==1)=nan;
        %combine with what we already have
        nans_mask2=nans_mask2.*nans_mask_swap;
        %exclude more?
        prompt = 'Select regions to exclude? y/n ';
        exclude_flag = input(prompt,'s');
    end
    
    %Center Selection, next select the center:
    disp('Select pattern center')
    figure(2)
    [x_center,y_center]=ginput(1)
    
end

%In not user input centers, look them up:
if center_lookup_on==1
    x_center=double(mie_x0);
    y_center=double(mie_y0);
end

x_0=int16(x_center);
y_0=int16(y_center);%making these unsigned, screws up the conditions down stream,
%need to be signed to allow for negative values when subtracted. 
%toc %c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Angle Matrices: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Matrices of theta and phi angles for the pupil are calculated externally
%by generate_angle_matrices and read in

thetas_new=zeros(size(img));

%disp('e- revisit this.. going to depend on index of media no?')
%tic

%we can correct for angle compression here, replacing a previous more 
%complicated algebraic method:
%happens here as it requires thetas. 
%hmm this is going to be different with water vs glycerine for example no? 
angle_compression_matrix=cosd(asind(1.33/1.515*sind(thetas)));
angle_compression_matrix(isnan(angle_compression_matrix))=1;
img=img.*angle_compression_matrix;
%toc

%plot angle space here:
%{
tic
xs=thetas.*cos(phis);
ys=thetas.*sin(phis);
z=img;
figure(13)
surf(real(xs),real(ys),real(z),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
view(0,90)
axis image
toc
disp('plot flattened angle space')
%}

%disp('g - perhaps remove this functionality..')
%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's loop over shifted transform points:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is almost legacy, may take this functionality out...

counter=1;
if do_looping==1
    x_extreme_neg=-20;
    x_extreme_pos=20;
    y_extreme_neg=-20;
    y_extreme_pos=20;
else
    x_extreme_neg=0;
    x_extreme_pos=0;
    y_extreme_neg=0;
    y_extreme_pos=0;
end    
    
%for runs=1:2
runs=1;
while runs<=2    
    if hard_coded_on==1
        runs=2;
        thumbs_up_down='1';
    end
for x_shift=x_extreme_neg:5:x_extreme_pos    
    for y_shift=y_extreme_neg:5:y_extreme_pos    
        
     x_0_mod=x_0+x_shift;
     y_0_mod=y_0+y_shift;

%toc %g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('h - this is the angle matrices for the off-center case, and the transform, seems pretty fast considering')
%tic

%first, find the angles for our new origin point, phi_prime and theta_prime:
y_fov2=double(y_0_mod-cir_y);%so this is the x and y location of pattern center...
x_fov2=double(x_0_mod-cir_x);
argument2=sqrt(x_fov2^2+y_fov2^2);%so this is the r,
%if argument2<cir_rad %w/ old cir radius
if argument2<cir_rad_to_mask %w/ newer smaller circle mask, if it falls inside mask circle
    %interp method finds theta of that center point. 
    theta_prime=interp1(pp_coord,angle_out,argument2/(cir_rad*1.0),'pchip');
    %need to wrap?..for some reason, this doesn't work:    
    %phi_prime=atan2(y_fov,x_fov);            
    if x_fov2>0
        phi_prime=atan(y_fov2/x_fov2);
    else
        phi_prime=atan(y_fov2/x_fov2)+pi;
    end
else %if it falls outside the circle, get nans:
    theta_prime=nan;
    phi_prime=nan;
end

theta_prime;%Note: if you want the patter center theta angle out for debugging, here it is...
phi_prime;

%{
%indicated the center.
figure(13)
hold on
plot3(theta_prime*cos(phi_prime),theta_prime*sin(phi_prime),1,'rX')
hold offtoc
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's explore rotation matricies, this gives the correct transmormation
%after previously working with a misleading (and wrong) just subtraction of
%vectors.

%first we need x, y, and z for the thetas...  
x_matrix=cos(phis).*sind(thetas);
y_matrix=sin(phis).*sind(thetas);
z_matrix=cosd(thetas);

%diagnostic figures of angle matrices:
%{
figure(20)
imagesc(thetas)
figure(21)
imagesc(phis)
figure(22)
imagesc(x_matrix)
figure(23)
imagesc(y_matrix)
figure(24)
imagesc(z_matrix)
%}

%next transform with rotation matrices, the below is after some trial and
%error, some confusion about which axis is which. Math from polar trig..
beta=atand(tand(theta_prime)*cosd(phi_prime*180/pi));
alpha=acosd(cosd(theta_prime)/cosd(beta));

if phi_prime>0 && phi_prime<pi
    alpha=-alpha;
end

%gamma=45;%z-rotation
gamma=0;

Rx=[1 0 0;0 cosd(alpha) sind(alpha);0 -sind(alpha) cosd(alpha)];
Ry=[cosd(beta) 0 -sind(beta);0 1 0;sind(beta) 0 cosd(beta)];
Rz=[cosd(gamma) sind(gamma) 0;-sind(gamma) cosd(gamma) 0;0 0 1];

%R=Rx*Ry*Rz;
R=Rz*Rx*Ry;
%R=Rz*Ry*Rx;

x_prime_matrix=R(1,1).*x_matrix+R(1,2).*y_matrix+R(1,3).*z_matrix;
y_prime_matrix=R(2,1).*x_matrix+R(2,2).*y_matrix+R(2,3).*z_matrix;
z_prime_matrix=R(3,1).*x_matrix+R(3,2).*y_matrix+R(3,3).*z_matrix;

%some diagnostic images of transformed versions:
%{
figure(32)
imagesc(x_prime_matrix)
figure(33)
imagesc(y_prime_matrix)
figure(34)
imagesc(z_prime_matrix)
%}

%transformed thetas
thetas_new=atan2d(sqrt(x_prime_matrix.^2+y_prime_matrix.^2),z_prime_matrix);

%need new phis for the mapping as well:
phis_new=atan2(y_prime_matrix,x_prime_matrix);
%toc h

%disp('i_1')
%tic
%figure(195)
%imagesc(phis_new)
%figure(196)
%imagesc(thetas_new)

%OK> so here is where we have access to the transformed thetas and phis..
%if we spit out these matrices, we could stitch together angular maps and
%then analyze in angle space...? 

%trying to visualize.... in 3D the scattering phase function. 
%{
Iblur = imgaussfilt(img, 10);
%Iblur = imgaussfilt(img, 15);
%[x_todisp,y_todisp,z_todisp]=sph2cart(phis_new,thetas_new*pi/180-pi/2,log10(Iblur+1));
[x_todisp,y_todisp,z_todisp]=sph2cart(phis_new,thetas_new*pi/180-pi/2,(Iblur+1));

%downsample so you can actually look at it: 
dn_sample_no=8;
x_todisp_ds=x_todisp(1:dn_sample_no:end,1:dn_sample_no:end);
y_todisp_ds=y_todisp(1:dn_sample_no:end,1:dn_sample_no:end);
z_todisp_ds=z_todisp(1:dn_sample_no:end,1:dn_sample_no:end);
r_todisp=log10(Iblur+1);
r_todisp=ones(size(Iblur));
r_todisp_ds=r_todisp(1:dn_sample_no:end,1:dn_sample_no:end);
%}

%disp('hey debugging here:')
%disp('max/min of phi')
%max(max(phis_new))
%min(min(phis_new))
%disp('max/min of theta')
%max(max(thetas_new_2))
%min(min(thetas_new_2))
%disp('max/min of img')
%max(max(img))
%min(min(img))
%phi_prime
%theta_prime
%alpha
%beta
%disp('end of debugging')

%{
%clf(19) %3D visualization
figure(23)
%figure
hold on
surf(x_todisp_ds,y_todisp_ds,z_todisp_ds,r_todisp_ds);%,'linestyle','none');
%surf(x_todisp,y_todisp,z_todisp);
az = 90;
el = -90;
colormap('autumn')
view(az, el);
axis image
%}

%so we either look at only the pixels selected, defined by nans_mask2:
    
if hard_coded_on~=1
    thetas_new=thetas_new.*nans_mask2;
    phis_new=phis_new.*nans_mask2;%added phi functionality
    mask_out=nans_mask2;
    sat_mask=sat_mask.*nans_mask2;
else
    %or by using external mask, if we want ti the same every time
    thetas_new=thetas_new.*mask_in;
    phis_new=phis_new.*mask_in;%added phi functionality
    mask_out=mask_in;
    sat_mask=sat_mask.*mask_in;
    %mask_out=0;
end

%figure(69)
%imagesc(sat_mask)
%drawnow

%reshape images as lists: 
thetas_prime_list=reshape(thetas_new,1,[])';
phis_prime_list=reshape(phis_new,1,[])';%added phi functionality
values_list=reshape(img,1,[])';
sat_list=reshape(sat_mask,1,[])';

%hmm export an object that has the same info but not separate imgs
lists_out(:,1)=values_list;
lists_out(:,2)=sat_list;
lists_out(:,3)=thetas_prime_list;
lists_out(:,4)=phis_prime_list;

%figure(20)
%plot(thetas_prime_list,values_list)

%max(max(thetas_prime_list));
%min(min(thetas_prime_list));
%max(max(phis_prime_list));
%min(min(phis_prime_list));
%toc %i_1

%disp('i_2')
%tic
%discard nans, i.e. stuff that is outside the FOV

if speed_mode~=1

isnan_list=isnan(thetas_prime_list);
values_list(isnan_list)=[];
phis_prime_list(isnan_list)=[];
sat_list(isnan_list)=[];
thetas_prime_list(isnan_list)=[];
%toc %i_2

%figure(19)
%histogram(smoothed_sorted_values)
%figure(18)
%semilogy(sorted_thetas_list,smoothed_sorted_values)

%binning makes a downsampled version, this has become the default:
    [theta_vector,phi_vector, radial_avg_downsampled_matrix]=bin_it3c(thetas_prime_list,phis_prime_list,values_list,sat_list,how_many_phi_bins);
    radial_avg=(radial_avg_downsampled_matrix(end,:));
    if how_many_phi_bins~=1
    %does this still work even?..     
       beta_vector=bin_it_intensity(thetas_prime_list,phis_prime_list,values_list,sat_list);
    end
else
    theta_vector=0;
    radial_avg_downsampled_matrix=0;
      
end%speed mode    
    
    
%disp('k')
%tic

%defunct
%{
%g integrals: This is defunct code, instead of calculating gs in mie
%analysis routine, we do this externally now in conjunction with fits.
%Doesn't really make sense to calculate the integrals for a PARTIAL trace
if how_many_phi_bins~=1
gs = g_integral( theta_vector,phi_vector,radial_avg_downsampled_matrix );
%figure(78)
%plot(phi_vector,gs);
g_avg=mean(gs);
%title(num2str(g_avg))
%figure(24)
%plot(gs)
%hold on
else
    g_avg=0;
end
%}    

%if runs==1%&&do_looping==1
 %disable this
if 1==0 
    figure(67)
    %plot(theta_vector,radial_avg*1/8,'r')
    plot(theta_vector,radial_avg,'r')
    drawnow
    hold on

    disp('If it looks good press 1, otherwise press anything else: ')
    w = waitforbuttonpress;
    if w
        p = get(gcf, 'CurrentCharacter');
        disp(p) %displays the character that was pressed
        %disp(double(p))  %displays the ascii value of the character that was pressed
    end    
    thumbs_up_down=strcmp('1',p);
    figure(67)
    plot(theta_vector,radial_avg,'b')
    drawnow
    hold on
end% if statement: depress if it looks good..
    
    thumbs_up_down=1;
    center_vicinity(counter,1)=x_shift;
    center_vicinity(counter,2)=y_shift;
    center_vicinity(counter,3)=thumbs_up_down;
    counter=counter+1;
    
    end%x_shift
end%y_shifts

runs=runs+1;

%center finding logic.. this is crude, consider stripping out..
%disp('end of trial part')
x_coords=center_vicinity(:,1);
y_coords=center_vicinity(:,2);
which_coords=center_vicinity(:,3);
x_coords=x_coords(which_coords>0);
y_coords=y_coords(which_coords>0);
    
%mie offsets:
x_offset=int16(mean(x_coords));
y_offset=int16(mean(y_coords));

%figure(17)
%plot(x_coords,y_coords,'bX')

x_extreme_neg=x_offset;
x_extreme_pos=x_offset;
y_extreme_neg=y_offset;
y_extreme_pos=y_offset;

runs;
end %trial and final runs
%end looping over shifted centers
%toc %k

%disp('l')
%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Plots/Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for l=1:100
    theta=l/100*2*pi;
    x(l)=10*cos(theta);
    y(l)=10*sin(theta);
    z(l)=1020;    
end
%}

%{
figure(16)
%ok let's plot theta space..
%surf(real(argx),real(argy),log10(real(img)),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
surf((argx),(argy),(log10(img)),'EdgeColor','none','LineStyle','none','FaceLighting','phong')

hold on

%surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    angle_limit=115;
    axis([-angle_limit angle_limit -angle_limit angle_limit 0 1024])
%caxis([0 .001])
view(-90,90)
%axis image
axis equal
disp('plot flattened angle space, re-centered')
%}

%{
figure(15)
%and also before fixing..
hold on
%imagesc((flipud(log10(img))));
imagesc(log10(img));
%caxis([0 .001])
set(gca,'YDir','reverse') %seem to have to do this to get proper display
circle(cir_x,1024-cir_y,cir_rad)
% mask circle:
%}


%plot some features on top of our raw image..
if figs_off~=1
    
    figure(2)
    plot(x_0_mod,y_0_mod,'yX')
    plot(x_0,y_0,'rX','MarkerSize', 12,'LineWidth',2)
    circle(cir_x,cir_y,cir_rad_NA)
    plot(cir_x,cir_y,'bX')
    drawnow
    
    %for ISO's
    halfwidth=.15;
    %ok this should be angles for 1.54um spheres forward..
    if fw_on==1
        if color==525
            angles=[16.7;32.5;47.3;62.6;77.6;94;110.4;127.4;127.4;127.4];
        else%red
            angles=[22.7;41.1;59.7;78.6;98.7;126.8;151.8;172.6;172.6;172.6];
        end
    else
        %ok this should be the angles for 1.54 spheres backward:
        if color==525
            angles=[5.8;23.3;38.7;52.9;69.7;86;102.1;117.6;117.6;117.6];
        else
            angles=[7.4;28.2;53.2;81.3;101.4;120.3;138.9;157.3;157.3;157.3];
        end
    end
    
    [row,col] = find((thetas_new>(angles(1)-halfwidth )& thetas_new<(angles(1)+halfwidth))|...
        (thetas_new>(angles(2)-halfwidth )& thetas_new<(angles(2)+halfwidth))|...
        (thetas_new>(angles(3)-halfwidth )& thetas_new<(angles(3)+halfwidth))|...
        (thetas_new>(angles(4)-halfwidth )& thetas_new<(angles(4)+halfwidth))|...
        (thetas_new>(angles(5)-halfwidth )& thetas_new<(angles(5)+halfwidth))|...
        (thetas_new>(angles(6)-halfwidth )& thetas_new<(angles(6)+halfwidth))|...
        (thetas_new>(angles(7)-halfwidth )& thetas_new<(angles(7)+halfwidth))|...
        (thetas_new>(angles(8)-halfwidth )& thetas_new<(angles(8)+halfwidth))|...
        (thetas_new>(angles(9)-halfwidth )& thetas_new<(angles(9)+halfwidth))|...
        (thetas_new>(angles(10)-halfwidth )& thetas_new<(angles(10)+halfwidth)));
    
    figure(2)
    %plot ISOs on top
    plot(col,row,'color',[.7 .7 .7],'marker','.','LineStyle','none')
end

%output mie center
mie_center_x=double(x_center+x_offset);
mie_center_y=double(y_center+y_offset);

mie_center(1,1)=mie_center_x;
mie_center(1,2)=mie_center_y;

%plot the theta map for debugging purposes..
%{
%figure(14), imagesc(thetas_new)
figure(14),imagesc(image_scaling_mask)
set(gca,'YDir','reverse')
colorbar
hold on
plot(xindex,yindex,'rX')
plot(cir_x,cir_y,'bX')
hold off
%}
%toc %-l
end %function end

%Z(isnan(Z))=0;
%TheMatrix(isnan(TheMatrix)) = [];

