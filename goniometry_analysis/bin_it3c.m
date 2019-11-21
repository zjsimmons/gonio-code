function [bins,phi_bins,bin_vals] = bin_it3c( angle_vector,phi_vector,val_vector,sat_vector,how_many_phi_bins )
%[theta_vector_downsampled,phi_vector_downsampled, radial_avg_downsampled]=bin_it3(theta_vector,phi_to_spit_out,radial_avg);
%DESCRIPTION: 
%function to bin/average traces as a fn of angle, this way raw traces that
%have different numbers of points are all the same size for use later,

%This program has grown from that original purpose to become a lot more
%complicated/include a lot more functionality. In addition to the binned,
%smoothed trace, scattering as fn of theta, it also outputs numerically
%integrated g (up to what NA you can see) as a fn of phi. This starts to
%get at information in the anisotropic tissue case. 

%General Workflow: 
%The image pixel value list (val_vector) and theta list (angle_vector) are
%sorted according to theta. Looping over a user-defined no of phi bins (pie
%pieces) in the BFP image, each sub-list is sorted according to theta and
%the trace generated for each slice. g (up to that limited NA) is computed
%for each slice as well. 

%There are a couple important user-defined parameters:
%sat_frac_threshold is the amount of a theta bin that can be saturated 
%before you decide to throw it out.  
%count_threshold - fewer than this no of theta values means you don't
%compute an average. 

%INPUTS: 
%val_vector - sorted list of pixels to include in trace
%angle_vector - list of corresponding thetas
%phi_vector - list of corresponding phis
%sat_vector - list of corresponding saturated pixels

%OUTPUTS: 
%bins - theta bins
%phi_bins - phi bins
%bin_vals - average at that bin
%gs -this is a vector of g values numerically integrated from the trace


%FUNCTIONS CALLED: 

%NOTES:

%update 4-28-2016, added functionality to deal with saturation..

%updated 7-17-2016, added functionality to incorporate phi as well.. this
%allows for returning a numerically computed value of g as a function of
%phi as well, has some problems, namely you don't have the whole phase
%function, so it's not really accurate, but it can show the shape in g as
%you go around the phase function if the scattering is not axially
%symmetric, which is interesting.

%should I add functionality to smooth as well/fill in gaps?? - need to
%think about this, discard then smooth? how much is suitable amount of
%saturation..?

%zjs, cleaned up 10-27-2016

%updated 1-17-2018, changed main logic to use accumarray, much faster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('j_1')
%tic
%turn on g(phi) stuff?
g_phi_on=0;

%since this version looks at phi bins, we need to re-sort accordingly: 
[sorted_phis,index_list]=sort(phi_vector);
sorted_thetas_by_phi=angle_vector(index_list);
sorted_img_by_phi=val_vector(index_list);
sorted_sat_by_phi=sat_vector(index_list);

%re: theta
angle_step_size=0.1;
max_angle=180;
bins=zeros(1,180/angle_step_size+1);
bin_vals=zeros(how_many_phi_bins,180/angle_step_size+1);

%re: phi:
%how_many_phi_bins=16;%choose how many phi bins (pie slices you want)
%how_many_phi_bins=32;
%how_many_phi_bins=1;
phi_step_size=2*pi/how_many_phi_bins;
%max_phi=359;

%max_pixel_val=1023;
%wait where does this no come from?- would be good to pin this down, rather
%than just rely on it being the brightest pixel., i.e. is there saturation?
%max_pixel_val=max(max(val_vector))
%min_pixel_val=min(min(val_vector));

sat_frac_threshold=0.01;% this is the amount of a theta bin that can be
%saturated before you decide to throw it out. Would be nice to have a
    %better method for choosing this number. 1% is a round no. 
sat_frac_threshold=1;%Note: to turn the sat-thresholding off, set to 1, NOT 0 
count_threshold=10;% this is how many counts you need in a theta bin to 
    %calculate an average
%count_threshold=0;
    
%toc %j_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Main Logic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sat_bin_count=0;
%loop over phi bins
phi_counter=1;
%for phi_bin=(0:phi_step_size:max_phi)*2*pi/180
for phi_bin=(-pi:phi_step_size:pi-.001)%looping over phi.. 
    %disp('j_2')
    %tic    
    phi_bins(phi_counter)=phi_bin;%so this is the phi bin value at hand
    
    %select out slice of phi: thetas, phis, pixel vals
    pie_slice_theta=sorted_thetas_by_phi(sorted_phis>=phi_bin & sorted_phis<(phi_bin+phi_step_size));
   % pie_slice_phi=sorted_phis(sorted_phis>=phi_bin & sorted_phis<(phi_bin+phi_step_size));
    pie_slice_img=sorted_img_by_phi(sorted_phis>=phi_bin & sorted_phis<(phi_bin+phi_step_size));
    pie_slice_sat=sorted_sat_by_phi(sorted_phis>=phi_bin & sorted_phis<(phi_bin+phi_step_size));
    
    %how many pixels per pie slice...
    %phi_bin_length=length(pie_slice_theta)
    
    %this lets you see the 'pie slice'
    % xs=sin(pie_slice_phi).*pie_slice_theta;
    % ys=cos(pie_slice_phi).*pie_slice_theta;
    % pause(.4)
    % figure(90)
    % hold on
    % scatter(xs,ys)
       
    counter=1; %theta count IN a phi slice
       
    %Note: want to limit the total theta so as to make them all the
    %same, perhaps if different pie slices go to different thetas, it's
    %skewing the result- this is a problem when off-axis...  
    %toc -j_2
    
    %disp('j_3')
    %tic
    %so this is the slow bit, now done with accumarray
    %so for this we need subscripts rather than a range, note: subs can't be zero
    
    pie_slice_theta2=uint16(floor(pie_slice_theta*10))+1;
    sz=[1801 1];
    
    %so this finds the mean of img values all subscripted by the same theta-related values.. 
    %subscripts are in pie_slice_theta2, vals in pie_slice_img, sz forces
    %the output size to be the full angle-range, @mean tells it to
    %calculate the mean of those subsets, nan means entries with no values
    %are padded with nan    
    
    %main operation: 
    bin_vals(phi_counter,:) = accumarray(pie_slice_theta2,pie_slice_img,sz,@mean,nan);    
    bins=0:.1:180;
    
    %Some logic to make the result more robust: 

    %1) logic to exclude if count is too small...
    %also want the count.. don't want to include elements that don't have
    %many to average.. but this probably isn't really a problem.. only
    %going to happen at the very center and at the edge..
    %bin_counts(phi_counter,:) = accumarray(pie_slice_theta2,pie_slice_img,sz,@sum,nan);    
    %
    bin_counts(phi_counter,:) = accumarray(pie_slice_theta2,1,sz);    
    %figure(100)
    %plot(bin_counts)
    bin_vals(bin_counts<count_threshold)=nan;
    %}
    
    %2) functionality to exclude as a function of theta based on fraction
    %saturated: 
    pie_slice_sat(isnan(pie_slice_sat))=0;
    bin_sat_frac(phi_counter,:) = accumarray(pie_slice_theta2,pie_slice_sat,sz,@mean,nan);
    %or one could do an absolute count of number of saturated.. 
    %bin_sat_frac(phi_counter,:) = accumarray(pie_slice_theta2,pie_slice_sat,sz,@sum,nan);
   
    sat_frac_threshold=1;
    bin_vals(bin_sat_frac>sat_frac_threshold)=nan;  
   
   %{
    figure(101)
    plot(bin_sat_frac)
    title('bin sat frac')
        
   % bin_val_size=size(bin_vals)
    figure(102)
    plot(bin_vals)
    title('bin vals')
    
    figure(103)
    plot(bin_counts)
    title('bin counts')
    %}
    
    
    %toc %j_3
    
%}
%for debugging, spits out the number of saturated bins:
%disp(strcat('there were_',' ',num2str(sat_bin_count),' saturated bins'))
%ok so we have bin_vals(n,m) where n is the phi angle, m is the theta angle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Re: g(phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Next let's extract some phi-dependent info, this needs improvement.
%{
if g_phi_on==1
    % first we need to normalize each. Q: what about the saturated fw
    % scattering and the fact we don't have the bw piece?
    
    %whatisthemin=min(min(bin_vals))% zero evidently
    %isthereanynans=isnan(bin_vals)
    bin_vals(isnan(bin_vals))=0;
    %above is the least sophisticated way of dealing with missing spots, will
    %probably need to be updated.
    
    %let's make the large angles not zero and see what happens to the g's it
    %spits out:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% improve this? but how to choose a good background to add?
    %%%%%%%%%%%%%%%%%%%%%%%
    theta_padding=ones(how_many_phi_bins,counter-1)*.1; %hmm this is very crude, must be a better way, don't want to add to the whole thing...
    bin_vals=bin_vals+theta_padding;
    
    %actually calculating g:
    % 1) normalization of p(\theta) in phi bins
    %note, transpose each so columns are theta
    integrands = bsxfun(@times,bin_vals',sind(bins)')*2*pi*angle_step_size;
    integrals=sum(integrands);
    normalization_factors=1./integrals; %1xno bins row % norm factors for each phi bin
    bin_vals_normalized=bsxfun(@times,bin_vals,normalization_factors');%this multiplies
    %each all the rows - p(thetas) by the appropriate norm factor for that
    %row (phi bin)
    
    %next calculate g's for each phi bin, multiply the p(thetas) by the
    %appropriate stuff to get the integral for g's:
    integrands2 = bsxfun(@times,bin_vals_normalized',(sind(bins).*cosd(bins))')*2*pi*angle_step_size;
    gs=sum(integrands2);
    
    %to display, wrap around
    padded_gs=[gs gs(1)];
    phis_to_display=0:360/how_many_phi_bins:360;
    
    figure(31)
    %subplot(2,1,1)
    hold on
    plot(phis_to_display,padded_gs)
    xlim([0 360])
    ylim([.4 1])
    xlabel('\phi')
    ylabel('g')
    title('g(\phi) computation brain matter \phi dependent g')
    maxg=max(gs)
    ming=min(gs)
    spread=(ming-maxg)/mean(gs)
    
%figure(15)
%polar(phis_to_display*pi/180,padded_gs)

%figure(13)
%plot(bins, bin_vals)
%plot(bins, integrands)
%plot(bins, bin_vals_normalized)

else
    gs=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  End of Re: g(phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

end

