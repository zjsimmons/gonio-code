function fw_portion=stitcher3_fn(lists_out,bw_flag)

%DESCRIPTION: This is a minor routine for stacking oriented acquisitions, 
%returning a binned version to facilitate averaging. 

%so this script will grow into a routine to stack/stitch angle-flattened
%image output. i.e. make actual 4pi 


%INPUTS: 
    %lists_out - this is the grid from the main analysis program that has:
    %brightness|saturation|thetas|phis
    %bw_flag - tells if this is fw or bw facing. 

%OUTPUTS: 
    %fw_portion is the binned (in phi and theta) output

%FUNCTIONS CALLED: 
    %mie_calc_fn.m (and a_n_b_n_fn.m and pi_tau_fn.m) - mie plotting
    %mmpolar.m - function to plot coverage plots

%zjs, cleaned up 8-15-2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       binned grid generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a couple of additional robust-ness checks: 
sat_threshold=0.01;
%sat_threshold=1;%Note: to turn the sat-thresholding off, set to 1, NOT 0
count_threshold=50;


%note: normalization is arbitrary, shouldn't do anything to g calculations
%for example, but it is nice to match it to Mie for example, to give a
%visual comparison. This is for getting sucha  factor: 
norm_factor=1;
%norm_factor=mie_calc_fn(.63,1.54)%norm factor is to normalize to phase 

%bw_flag=0;
fwqmax=7;
bwqmax=4;
fw_scaling_factor=norm_factor/2;
bw_scaling_factor=norm_factor/20;
fw_scaling_factor=1;
bw_scaling_factor=1;

%create angle grid: 
%phi x theta, 
angle_res=1;
no_phi=360/angle_res; no_theta=180/angle_res;

%ok so we have a big grid where the columns are img/saturated?/thetas/phis
if bw_flag==1
   qmax=bwqmax;
else
   qmax=fwqmax;
end

for q=1:qmax
angle_grid=zeros(no_phi,no_theta);
count_grid=angle_grid;
sat_count_grid=angle_grid;
list_at_hand=lists_out{q};

if bw_flag==1
    scaling_factor=bw_scaling_factor;
    thetas_to_examine=180-list_at_hand(:,3);     
else
    scaling_factor=fw_scaling_factor;
    thetas_to_examine=list_at_hand(:,3);
end

img_to_examine=list_at_hand(:,1)*scaling_factor; 
phis_to_examine=list_at_hand(:,4);
sats_to_examine=list_at_hand(:,2);
sats_to_examine_save=sats_to_examine;

%kill all saturated and nan values first? 
%{
sats_to_examine(sats_to_examine==1)=nan;
sats_to_examine(sats_to_examine==0)=1;

%ditch all the excluded, so we don't have to sort through them: 
%multiply by nans mask.. 
%note this will be changed to instead exclude by a fraction at a particular
%theta that is saturated, rather than just killing all saturated pixels. 
%
img_to_examine=img_to_examine.*sats_to_examine;
phis_to_examine=phis_to_examine.*sats_to_examine;
thetas_to_examine=thetas_to_examine.*sats_to_examine;
%thetas_to_examine_save=phis_to_examine;
%}

%discard nans- this excludes pixels outside the NA circle.. 
sats_to_examine_save=sats_to_examine_save(~isnan(phis_to_examine));
%img_to_examine=img_to_examine(~isnan(img_to_examine));
img_to_examine=img_to_examine(~isnan(phis_to_examine));
phis_to_examine=phis_to_examine(~isnan(phis_to_examine));
thetas_to_examine=thetas_to_examine(~isnan(thetas_to_examine));

%hmm selecting the groups at angles is very slow, look through the lists
%instead, sort to angles: 

for p=1:length(img_to_examine)

    val_at_hand=img_to_examine(p);
    theta_at_hand=thetas_to_examine(p);
    phi_at_hand=phis_to_examine(p);
    %pull out the locations:
    theta_at_hand_index=uint32(floor(theta_at_hand/angle_res))+1;
    phi_at_hand_index=uint32(floor((phi_at_hand+pi)*180/pi/angle_res))+1;
    %sort into the appropriate angle bins: 
    angle_grid(phi_at_hand_index,theta_at_hand_index)=angle_grid(phi_at_hand_index,theta_at_hand_index)+val_at_hand;
    count_grid(phi_at_hand_index,theta_at_hand_index)=count_grid(phi_at_hand_index,theta_at_hand_index)+1;
   
    %also want a sat count to find the fraction saturated.. 
    sat_at_hand=sats_to_examine_save(p);
    sat_count_grid(phi_at_hand_index,theta_at_hand_index)=sat_count_grid(phi_at_hand_index,theta_at_hand_index)+sat_at_hand;
        
end

%some debugging figures:
%
%count trace.. 
count_trace=sum(count_grid,1);
%figure(100)
%plot(count_trace)
%hold on
%sat trace
sat_trace=sum(sat_count_grid,1);
%figure(100)
%plot(sat_trace)
%drawnow
%hold off
%pause(.5)
%size(sat_trace)
%}

%let's institute a couple protections: 
%1) logic to throw out thetas that are saturated: 
exclude_bc_sat=sat_trace./count_trace;
%figure(101)
%plot(exclude_bc_sat)
%drawnow
%hold on
exclude_bc_sat(exclude_bc_sat>sat_threshold)=nan;
%plot(exclude_bc_sat)
%title('saturated fraction')
%drawnow
%hold off

%2) logic to throw out theta b/c of low counts: count_threshold
exclude_bc_counts=ones(size(count_trace));
%figure(102)
%plot(exclude_bc_counts)
%hold on
exclude_bc_counts(count_trace<count_threshold)=nan;
%plot(exclude_bc_counts)
%hold off
exclude_combined=ones(size(count_trace));
exclude_combined(isnan(exclude_bc_counts)|isnan(exclude_bc_sat))=nan;

%figure(103)
%plot(exclude_combined)

%sat_img=reshape(sats_to_examine_save,2048,[]);
%figure(102)
%imagesc(sat_img)
%drawnow

%theta_img=reshape(thetas_to_examine_save,2048,[]);
%figure(103)
%imagesc(theta_img)
%drawnow

%pause

%find avg- so this is the trace averaged over phi
avg_angle_grid=angle_grid./count_grid;

%size(exclude_combined)
%size(avg_angle_grid)

avg_angle_grid=avg_angle_grid.*exclude_combined;
%saving/compiling a stack for all the images (after binning/sorting)
if bw_flag==1
%    avg_angle_grids(:,:,fwqmax+q)=avg_angle_grid;%hmm..
    avg_angle_grids(:,:,q)=avg_angle_grid;
else
    avg_angle_grids(:,:,q)=avg_angle_grid;
end    

%{
avg_trace=nanmean(avg_angle_grid);
figure(123)
plot(avg_trace)
hold on
drawnow
%}

%ok cool, we should be able to next plot multiple binned slices
%corresponding to different angle ranges, then average them together..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               angular coverage figure: flat projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[theta_grid,phi_grid]=meshgrid(1:angle_res:180,1:angle_res:360);

theta_grid_list=reshape(theta_grid,1,[])';
phi_grid_list=reshape(phi_grid,1,[])';
r_grid_list=reshape(avg_angle_grid,1,[])';

theta_grid_list=theta_grid_list(~isnan(r_grid_list));
phi_grid_list=phi_grid_list(~isnan(r_grid_list));
%r_grid_list=reshape(avg_angle_grid,1,[])';


%angle coverage plot..
%{
figure(33)
colors={'.r','.b','.g','.c','.y','.k'};
colors=[colors colors];
color_selection=colors{q};

%perimeter
mmpolar((1:360)*pi/180,180*ones(1,360),'.b')
hold on
%coveraged
mmpolar(phi_grid_list*pi/180,theta_grid_list,color_selection)
title('4pi coverage plot')
drawnow
%}

end

fw_portion=avg_angle_grids(:,:,1:qmax);
fw_portion_flattened=nanmean(fw_portion,3);
flattened=fw_portion_flattened;
flattened_avg_trace=nanmean(flattened,1);
    %figure(18)
    %plot(flattened_avg_trace,'b')

end


%%
%old:
%ok so we have a datacube of phi,theta,scan angle, need to flatten over
%scan angle, so the combined:  
%{
fw_portion=avg_angle_grids(:,:,1:7);
fw_portion_flattened=nanmean(fw_portion,3);
bw_portion=avg_angle_grids(:,:,8:11);
bw_portion_flattened=nanmean(bw_portion,3);


flattened_avg_angle_grids=nanmean(avg_angle_grids,3);
[theta_grid,phi_grid]=meshgrid(1:angle_res:180,1:angle_res:360);

theta_grid_list=reshape(theta_grid,1,[])';
phi_grid_list=reshape(phi_grid,1,[])';
r_grid_list=reshape(flattened_avg_angle_grids,1,[])';

theta_grid_list=theta_grid_list(~isnan(r_grid_list));
phi_grid_list=phi_grid_list(~isnan(r_grid_list));


figure(34)
%polar(phi_grid_list*pi/180,theta_grid_list,color_selection)
mmpolar(phi_grid_list*pi/180,theta_grid_list,color_selection)

title('4pi coverage plot')
drawnow

%
%can also plot the averaged over phi trace: 
flattened_avg_trace=nanmean(flattened_avg_angle_grids,1);
figure(18)
plot(flattened_avg_trace*.75)

%%

%want a 3d plot to show the dimples for the dipole w/ linear polarization: 

figure(40)
polarplot3d(30*log10(flattened_avg_angle_grids'),'AngularRange',[0 2*pi],'RadialRange',[0 180],'PolarGrid',{6 8},'TickSpacing',0);
%polarplot3d(log10(flattened_avg_angle_grids'),'AngularRange',[0 2*pi],'RadialRange',[0 180]);

%polarplot3d((fw_portion_flattened'),'AngularRange',[0 2*pi],'RadialRange',[0 180]);
hold on
%caxis([-3 -1])
az = -180;
el = 90;
view(az, el);
axis image
set(gca, 'XTick', []);
set(gca, 'YTick', []);

%circle_black(0,0,180)
%%
polarplot3d(log10(avg_angle_grids(:,:,3)'),'AngularRange',[0 2*pi],'RadialRange',[0 180],'PolarGrid',{6 8},'TickSpacing',0);
%hold on
%polarplot3d(log10(avg_angle_grids(:,:,8)'),'AngularRange',[0 2*pi],'RadialRange',[0 180],'PolarGrid',{6 8},'TickSpacing',0);




%%
%r_grid_list(~isnan(r_grid_list))=1;
%{

[x_todisp,y_todisp,z_todisp]=sph2cart(phi_grid_list,theta_grid_list,r_grid_list);

figure(23)
hold on
surf(x_todisp,y_todisp,z_todisp,[0 1]);%,'linestyle','none');
%surf(x_todisp,y_todisp,z_todisp);
az = 90;
el = -90;
colormap('autumn')
view(az, el);
axis image

%}

%}