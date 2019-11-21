%
% This contains PS in H2O for calibration. I broke the script into 2, 1 for
% calibration, choosing masks, centers, etc. and saving them.. 1 for doing
% the analysis on samples of interest. This should hopefully make it easier
% for others to use. 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       General Set-Up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Have this set up with 2 'modes', setting the calibration or just looking
% at it. If the set_calibration flag is set, then it will let you choose
% masks and centers, save the result, etc. otherwise it just displays the
% current

set_calibration=1;
set_calibration=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        A. Forward Set-up:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%angle mapping for H2O
load('xout_yout.mat')
    
%existing lung masks & centers:
load('../../Data/Lung/2019_8_26/lung_fw_masks_centers.mat')

%H20 8-26 set
file=dir(strcat('../../Data/Lung/2019_8_26/47/*bgfw.mat'));
filename=strcat('../../Data/Lung/2019_8_26/47/',file.name);

load(filename);
bg_datacube=mean(datacube,2);

%fw_datacube=datacube;

%bg_datacube=0;

%this pulls the series settings from the acquisition tag
series_settings=tag{3,2}
[how_many columns]=size(series_settings);

%bg_of_interest=0;

%%%%%%%%%%%%%%%%%%%%% Analsysis Settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image and mask characteristics:
analysis_settings.image_size=[2048 2048];

analysis_settings.analysis_circle_x=1000;
analysis_settings.analysis_circle_y=1024;
analysis_settings.analysis_circle_NA=1180;

analysis_settings.mask_circle_x=analysis_settings.analysis_circle_x;
analysis_settings.mask_circle_y=analysis_settings.analysis_circle_y;
analysis_settings.mask_r=1050;
analysis_settings.mask_r=1100;

%analysis choices:
analysis_settings.subtract_bg_surf=1;
analysis_settings.fw_on=1;%1=fw, 0=bw
analysis_settings.color=600;
analysis_settings.use_input_mask_on=1;
analysis_settings.center_lookup_on=1;
analysis_settings.how_many_phi_bins=1;
analysis_settings.figures_off=0;
analysis_settings.speed_mode=0;

pp_coord=xout;
angle_out=yout;

%generate the angle matrices that will be used in every iteration, based on
%the pupil-plane mapping, xout and yout, i.e. pp_coordinate to angle map
[thetas,phis] = generate_angle_matrices(analysis_settings,xout,yout);

%ancillary analysis inputs: 
ancillary_inputs.pp_coord=xout;
ancillary_inputs.angle_out=yout;
ancillary_inputs.thetas=thetas;
ancillary_inputs.phis=phis;

if set_calibration  
    %temp dummy inputs:
    for n=1:7,masks_fw_in{n}=0;end
    mie_centers_fw_in=zeros(7,2);
    analysis_settings.use_input_mask_on=0;
    analysis_settings.center_lookup_on=0;
end

disp('forward set-up time:')
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         A. Backward Set-Up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%angle mapping, H2O
load('xout_yout.mat')
   
load('../../Data/Lung/mie_centers_masks_bw_in.mat')

%part of set
file=dir(strcat('../../Data/Lung/2019_8_26/47/*bgbw.mat'));
filename=strcat('../../Data/Lung/2019_8_26/47/',file.name);

load(filename);
bw_bg_datacube=datacube;

%this pulls the series settings from the acquisition tag
series_settings=tag{3,2}
[how_many columns]=size(series_settings);

%bg_of_interest=0;

%%%%%%%%%%%%%%%%%%%%% Analsysis Settings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image and mask characteristics:

analysis_settings_bw.image_size=[1920 1080];
analysis_settings_bw.analysis_circle_x=1135;
analysis_settings_bw.analysis_circle_y=550;
analysis_settings_bw.analysis_circle_NA=505;
analysis_settings_bw.mask_circle_x=analysis_settings_bw.analysis_circle_x;
analysis_settings_bw.mask_circle_y=analysis_settings_bw.analysis_circle_y;
analysis_settings_bw.mask_r=445;

%analysis choices:
analysis_settings_bw.subtract_bg_surf=0;
analysis_settings_bw.fw_on=0;%1=fw, 0=bw
analysis_settings_bw.color=600;
analysis_settings_bw.use_input_mask_on=1;
analysis_settings_bw.center_lookup_on=1;
analysis_settings_bw.how_many_phi_bins=1;
analysis_settings_bw.figures_off=0;
analysis_settings_bw.speed_mode=0;

pp_coord=xout;
angle_out=yout;

%generate the angle matrices that will be used in every iteration, based on
%the pupil-plane mapping, xout and yout, i.e. pp_coordinate to angle map
[thetas,phis] = generate_angle_matrices(analysis_settings_bw,xout,yout);

%ancillary analysis inputs: 
ancillary_inputs_bw.pp_coord=xout;
ancillary_inputs_bw.angle_out=yout;
ancillary_inputs_bw.thetas=thetas;
ancillary_inputs_bw.phis=phis;

disp('backward set-up time:')
toc

if set_calibration  
    %temp dummy inputs:
    for n=1:4,masks_bw_in{n}=0;end
    mie_centers_bw_in=zeros(4,2);
    analysis_settings_bw.use_input_mask_on=1;
    analysis_settings_bw.center_lookup_on=1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Loop over files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%more curated mouse lung sets:
%date_folder='2019_10_18/';
%tissue_types_cell_grid={'27_3';'27_4';'28_5';'28_6';'28_11';'28_13'}

%tissue_type='28_6';

%for ttn=1:1
%    tissue_type=tissue_types_cell_grid{ttn}

%next we want to sort through all the files we have.. 
how_many_files=1

thetas=0:179;
theta_min=15;
theta_max=180;

%load weights
    %load(strcat('../../Data/bnr/bnr_tissue_realign/processed/',tissue_type,'_trace_weights.mat'));
    %weights=mean(g_traces_weights,2);
weights=1;


for file_n=1:how_many_files
 
% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%fw

%september sphere set: 
file=dir(strcat('../../Data/Lung/spheres_3/*spheresfw.mat'));
filename=strcat('../../Data/Lung/spheres_3/',file.name);

%Human Lung, 9-25 sets.. 
%file=dir(strcat('../../Data/Lung/2019_9_25/',tissue_type,'/','*lung_',tissue_type,'_',num2str(file_n),'fw.mat'));
%filename=strcat('../../Data/Lung/2019_9_25/',tissue_type,'/',file.name);

%file=dir(strcat('../../Data/Lung/',date_folder,tissue_type,'/','*lung_',tissue_type,'_',num2str(file_n),'fw.mat'));
%filename=strcat('../../Data/Lung/',date_folder,tissue_type,'/',file.name);

%oct set
%file=dir(strcat('../../Data/Lung/',date_folder,tissue_type,'/','*',tissue_type,'_',num2str(file_n),'fw.mat'));
%filename=strcat('../../Data/Lung/',date_folder,tissue_type,'/',file.name);


bg_datacube=0;
load(filename);
fw_datacube=datacube-bg_datacube;
%fw_datacube=bg_datacube; 

%bw
%september sphere sets: 

file=dir(strcat('../../Data/Lung/spheres_3/*spheresbw.mat'));
filename=strcat('../../Data/Lung/spheres_3/',file.name);


%Human Lung sets
%file=dir(strcat('../../Data/Lung/2019_9_25/',tissue_type,'/','*lung_',tissue_type,'_',num2str(file_n),'bw.mat'));
%filename=strcat('../../Data/Lung/2019_9_25/',tissue_type,'/',file.name);

%file=dir(strcat('../../Data/Lung/',date_folder,tissue_type,'/','*lung_',tissue_type,'_',num2str(file_n),'bw.mat'));
%filename=strcat('../../Data/Lung/',date_folder,tissue_type,'/',file.name);

%oct
%file=dir(strcat('../../Data/Lung/',date_folder,tissue_type,'/','*',tissue_type,'_',num2str(file_n),'bw.mat'));
%filename=strcat('../../Data/Lung/',date_folder,tissue_type,'/',file.name);

%
bw_bg_datacube=0;
load(filename);
bw_datacube=datacube;
bw_datacube=datacube-bw_bg_datacube;
%bw_datacube=bw_bg_datacube;

%}

before=datetime;
%loop through our colors

% examine the individual colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for color_n=1:4
    
    %plot the Mie pattern:
    lambda=series_settings(color_n,1)/1000; 
    %color specific masks, etc. 
    %load(strcat('../../Data/bnr/sphere_calibration/masks_centers_fw_',num2str(1000*lambda),'.mat'))
    % load(strcat('../../Data/Lung/2019_5_20/masks_fw.mat'))
    
    %IDS color specific backward masks
    %  load(strcat('../../Data/bnr/sphere_calibration/masks_centers_bw_',num2str(1000*lambda),'.mat'))
    
    % usenormalization procedure based on Mie.?
    %norm_factors(color_n)=mie_calc_fn_wcv(lambda,1.49,0.10,'silica','PDMS','circ');
    %fused silica in PDMS for Bn(r) phantom
    %norm_factors(color_n)=mie_calc_fn_wcv(lambda,.96,0.075,'silica','PDMS','circ');
    %polystyrene in H2O calibration..
    norm_factors(color_n)=mie_calc_fn_wcv(lambda,1.54,0.05,'polystyrene','H2O','circ');
    
    angle_cube_fw=squeeze(fw_datacube(color_n,:,:,:));
    angle_cube_bw=squeeze(bw_datacube(color_n,:,:,:));
    
    %forward analysis
    for angle_n=1:1:7
        angle_n    ;
        bg_of_interest=0;           
        img_of_interest=squeeze(angle_cube_fw(angle_n,:,:));
        
        %main function
        figure(2);hold off;
        [lists_out_angles{angle_n},mie_centers_fw(angle_n,:),mask_out_cell_fw{angle_n},r_vector,radial_avg_matrix]=...
        mie_analysis_fn_universal8(img_of_interest,bg_of_interest,masks_fw_in{angle_n},mie_centers_fw_in(angle_n,1),mie_centers_fw_in(angle_n,2),ancillary_inputs,analysis_settings);    
        
        figure(2)
       % pause(1)
        close(2)
        %visualize phase function
        %phase_fn_vis_3d(lists_out_angles{angle_n},'raw');
        %plot the trace
        
    end %angle no

    fw_scan=stitcher3_fn(lists_out_angles,0);
    %}
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization stuff for fw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %let's look at the phi/theta grids, kind of interesting visualization: 
    %{
    figure(100)
    for n=1:5
       subplot(1,5,n)
       imagesc(log10(fw_scan(:,:,n)))
        axis image
    end
    %}
    
    %interpolated version for visualization purposes.. 
   fw_scan_avg=nanmean(fw_scan,3);        
   %fw_trace_avg=nanmean(fw_scan_avg,1);
   %figure(18),plot(fw_trace_avg*norm_factor)
   
   %figure(18),plot(r_vector,radial_avg_matrix)

    %visualize phase function
   % phase_fn_vis_3d(fw_portion_flattened,'binned');
    %}
    
    %optimize to visualize..
    %{
    fw_portion_flattened=fw_scan_avg;
    fw_portion_flattened2=nanmean(fw_portion_flattened,1);    
    figure(20)
    plot(fw_portion_flattened2)
    hold on
    vs=fw_portion_flattened2(7:end);
    xs=7:180;
    interped=interp1(xs,vs,1:180,'linear','extrap')
    figure(20)
    plot(interped)
    interped_matrix=repmat(interped,360,1);
    phase_fn_vis_3d(interped_matrix,'binned');
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %backward analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
%    close(2)
    for angle_n=1:1:4
        angle_n;
        bg_of_interest=0;           
        img_of_interest=squeeze(angle_cube_bw(angle_n,:,:));
        %max_pixel_val=max(max(img_of_interest))
           
        %main function
        [lists_out_angles{angle_n},mie_centers_bw(angle_n,:),mask_out_cell_bw{angle_n},r_vector,radial_avg_matrix]=...
        mie_analysis_fn_universal8(img_of_interest,bg_of_interest,masks_bw_in{angle_n},mie_centers_bw_in(angle_n,1),mie_centers_bw_in(angle_n,2),ancillary_inputs_bw,analysis_settings_bw);    
        
        %visualize phase function
        %phase_fn_vis_3d(lists_out_angles{angle_n},'raw');
        %plot the trace
        %figure(18),plot(r_vector,radial_avg_matrix)
        % pause
    end %angle no
    
    bw_scan=stitcher3_fn(lists_out_angles,1);
    
    
    %interpolated version for visualization purposes.. 
    %bw_portion_flattened=squeeze(nanmean(bw_scan,1));
    %figure(19)
    %plot(thetas,bw_portion_flattened(:,1:4))
    
    %fw-bw-absolute matching, set manually
    fw_factors(1,1)=.33;
    fw_factors(2,1)=.33;
    fw_factors(3,1)=.2;
    fw_factors(4,1)=.15; 
  
    bw_factors(1,1)=.002;
    bw_factors(2,1)=.004;
    bw_factors(3,1)=.0005;
    bw_factors(4,1)=.0004; 
       
    dc_bg=0;
        
    %this is back-scattering ratio?
   %debug% scattered_back_nos(color_n,file_n)=scatter_back(bw_scan*bw_factors(color_n,1),5);
    
    scan(:,:,1:7)=(fw_scan(:,:,1:7)*fw_factors(color_n)*norm_factors(color_n))-dc_bg;

   %hmm can kill bw portion.. or not use all of it if they're poor quality 
    scan(:,:,8:11)=bw_scan*bw_factors(color_n)*norm_factors(color_n);
   % scan(:,:,8)=bw_scan(:,:,3)*bw_factors(color_n)*norm_factors(color_n);

    %average the whole mess
    scan_avg=nanmean(scan,3);    
    scans_datacube(file_n,color_n,:,:)=scan_avg;
    
    %regarding the g extraction: 
    %
    %this is to find the min angle to incorporate in the fit:
    trace_avg=nanmean(scan_avg,1);
    figure(18),plot(thetas,trace_avg)
    %
    non_nan=1;
    while isnan(trace_avg(non_nan))
        non_nan=non_nan+1;
    end
    theta_min=thetas(non_nan);
    
    %fw_scans_colors{color_n}=stitcher3_fn(lists_out_angles,0);
    size(thetas);
    size(trace_avg);
    weights=ones(size(thetas))';
    
    %[g_int(color_n,file_n),g_fit(color_n,file_n),phasefn_normalized_backscatter(color_n,file_n),r_squared(color_n,file_n)]=fitting_g_vs_calculating_it_fn4(thetas',theta_min,theta_max,trace_avg',0);
    [g_int(color_n,1),g_fit(color_n,1),phasefn_normalized_backscatter(color_n,1),r_squared(color_n,1)]=fitting_g_vs_calculating_it_fn4(thetas',theta_min,theta_max,trace_avg',weights);
    title('HG fit')
    
    g_traces(:,1)=thetas';
    g_traces(:,color_n+1)=trace_avg';
    %}
    
end %color no

%data cubes.. 
g_traces_cube(:,:,file_n)=g_traces;

%color|g fit | g int| normalized backscatter | r_squared
g_chart=[.575; .650; .725; .8];
g_chart(:,2)=g_fit;
g_chart(:,3)=g_int;
g_chart(:,4)=phasefn_normalized_backscatter; %relevant to OCT comparison
g_chart(:,5)=r_squared;

g_chart_cube(:,:,file_n)=g_chart;

file_n

after=datetime;
set_time=after-before

figure(18)
title(num2str(file_n))

%file for each.. 
tag{1,1}='Notes/description:';
%tag{1,2}=strcat('looking at tissue, 50um thick, ',tissue_type,' ,fixed');
%tag{1,2}=strcat('looking at human lung tissue, 20um thick, ',tissue_type,' ,fixed');
tag{1,2}=strcat('looking at calibration, blah blah blah');

%save(strcat('../../Data/Lung/2019_5_21/',tissue_type,num2str(file_n),'b.mat'),'g_traces','g_chart','tag');
%save(strcat('../../Data/Lung/2019_8_26/',tissue_type,num2str(file_n),'.mat'),'g_traces','g_chart','tag');
%save(strcat('../../Data/Lung/2019_9_25/',tissue_type,num2str(file_n),'.mat'),'g_traces','g_chart','tag');
%save(strcat('../../Data/Lung/',date_folder,tissue_type,'_',num2str(file_n),'.mat'),'g_traces','g_chart','tag');

end %file no

r_squareds=squeeze(g_chart_cube(:,5,:));
rs_max=max(max(r_squareds));
rs_min=min(min(r_squareds));
disp(['r squareds range from ',num2str(rs_min,2),' to ',num2str(rs_max,2)])


%save summary from all 
%save(strcat('../../Data/Lung/2019_9_25/',tissue_type,'_summary_weighted_fits.mat'),'g_traces_cube','g_chart_cube','tag');
%save(strcat('../../Data/Lung/',date_folder,tissue_type,'_summary_weighted_fits.mat'),'g_traces_cube','g_chart_cube','tag');

%end


%% output and save the masks and centers.. 


%after first run, to set masks and centers: 
%mie_centers_fw_in=mie_centers_fw;
%masks_fw_in=mask_out_cell_fw;


%after first run, to set masks and centers: 
%mie_centers_bw_in=mie_centers_bw;
%masks_bw_in=mask_out_cell_bw;
