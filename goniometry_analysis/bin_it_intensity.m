function [beta_vector] = bin_it_intensity(sorted_thetas_list,sorted_phis_list,smoothed_sorted_values,sorted_sat_list)

%DESCRIPTION:
%This function pulls out iso's of intensity, circles in axially symmetric
%scattering to more oblong shapes in anisotropic tissue scattering. Those
%rings/ovals are then fit with the fitting function of Shuaib and Yao, see:
%Equi-intensity distribution of optical reflectance in a fibrous turbid medium
%APPLIED OPTICS / Vol. 49, No. 5 / 10 February 2010

%note: this only works well when you're on center, perhaps need to add
%logic to help with off-center. Diagnostic figure showing the fit is
%useful. Usually you can tell pretty quick if the fit to the oval looks
%reasonable, kind of neat to look at as well. 

%INPUTS:
%sorted_thetas_list - list of thetas from main analysis program
%sorted_phis_list - list of phis from main analysis program
%smoothed_sorted_values - list of scattering from main analysis program
%sorted_sat_list - where the saturated pixels are

%OUPUTS:
%beta_vector - list of betas for the chosen intensities (fractions of max)
%this are axis ratios, list tells you how it changes with brightness of
%scattering. 

%zjs, cleaned up 12-9-2016

%bin_it_intensity

%discard saturated pixels, sat mask=1 where saturated, so
sorted_thetas_list=sorted_thetas_list(sorted_sat_list~=1);
sorted_phis_list=sorted_phis_list(sorted_sat_list~=1);
smoothed_sorted_values=smoothed_sorted_values(sorted_sat_list~=1);
%discard any remaining nans?

%so first let's choose some intensities, we'll divide the whole range into
%some no of slices, discard zero and max slices, not so useful..
min_pixel_val=min(min(smoothed_sorted_values))
max_pixel_val=max(max(smoothed_sorted_values))
range_in_inputs=(max_pixel_val-min_pixel_val);
how_many_contours=4;
intensity_step=range_in_inputs/(how_many_contours+1);
intensities=(min_pixel_val+intensity_step):intensity_step:(max_pixel_val-intensity_step/2)

%so for generating the rings (ovals) how close to the chosen intensity
%should it be, some fraction of the step size?
intensity_spread=intensity_step/100;

color_string{1}='.r';
color_string{2}='.g';
color_string{3}='.b';
color_string{4}='.r';
color_string{5}='.g';

%next select out the ovals:
for n=1:length(intensities)
    %snap don't even need the actual intensities.
    %oval_intensity=smoothed_sorted_values(smoothed_sorted_values>(intensities(n)-intensity_spread)& smoothed_sorted_values<(intensities(n)+intensity_spread));
    %thetas and phis:
    oval_thetas=sorted_thetas_list(smoothed_sorted_values>(intensities(n)-intensity_spread)& smoothed_sorted_values<(intensities(n)+intensity_spread));
    oval_phis=sorted_phis_list(smoothed_sorted_values>(intensities(n)-intensity_spread)& smoothed_sorted_values<(intensities(n)+intensity_spread));
    %sort these according to phi?
    [oval_phis_sorted,index_list]=sort(oval_phis);
    oval_thetas_sorted=oval_thetas(index_list);
    %next resample to get equally spaced
    oval_phis_interped=-pi:pi/1000:pi;
    try
        oval_thetas_interped=interp1(oval_phis_sorted,oval_thetas_sorted,oval_phis_interped,'linear','extrap');
    catch
        disp('problem with bin_it_intensity.m')
    end
    %figure(77),plot(oval_thetas_interped),drawnow
    %next we need to fit..
    %pf=fittype('1/(4*pi)*(1-g^2)/(1+g^2-2*g*cos(x*pi/180))^(3/2)*a+b',...
    % 'dependent',{'y'},'independent',{'x'},'coefficients',{'g','a','b'})
    theta=oval_thetas_interped';
    phi=oval_phis_interped';
    %hmm some thetas are nans?
    
    pf=fittype('((abs(cos(phi-phi_0))/theta_x)^q+(abs(sin(phi-phi_0))/theta_y)^q)^(-1/q)',...
        'dependent',{'theta'},'independent',{'phi'},'coefficients',{'phi_0','theta_x','theta_y','q'});
    
    options = fitoptions(pf);
    options.Lower = [-pi 0 0 1];
    options.Upper = [pi 100 100 2];
    
    [myfit,gof]=fit(phi,theta,pf,options)
    
    %there must be an easier way to extract this stuff:
    phi_0=myfit.phi_0;
    theta_x=myfit.theta_x;
    theta_y=myfit.theta_y;
    q=myfit.q;
    oval_fit = @(phi) ((abs(cos(phi-phi_0))./theta_x).^q+(abs(sin(phi-phi_0))./theta_y).^q).^(-1./q);
    
    oval_fit_vector=oval_fit(oval_phis_interped);
    
    r_squared=gof.rsquare;
    
    %figure(71)
    %plot(myfit,phi,theta)
    %hold on
      
    %ok, for now let's have it polar plot
    %polar(THETA,RHO)
    figure(74)
    polar(oval_phis_interped,oval_thetas_interped,color_string{n})
    hold on
    polar(oval_phis_interped,oval_fit_vector,'k')
    title('scattering constant intensity isos for beta calc')
    
    beta=theta_y/theta_x;
    
    %x and y can be swapped so:
    if beta<1
        beta=1/beta;        
    end
    beta_vector(n,1)=beta;
    r_squared_vector(n,1)=r_squared;      
end

beta_vector
%figure(75)
%plot(beta_vector)
%ylim([0 3])
r_squared_vector

end

