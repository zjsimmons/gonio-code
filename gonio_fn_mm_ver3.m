%function []=gonio_fn_mm_ver2(savename,bw_flag)
function []=gonio_fn_mm_ver3(savename,mode_flag)

%% Program for scanning galvos, laser color, and capturing images for the
%goniometer experiment

% Requires data acquisition toolbox for NIDAQ control
% Requires micromanager and associated adapters for selected cameras

% This program communicates with the galvos via NIDAQ control toolbox,
%the laser via telegram messages, see: laser_fns2.m, and the cameras via
%micromanager, see camera_handler.m. Low-level logic is contained in the
%above .m fns (which return function handles to laser and camera functions)
%to maintain this as a simpler high-level program, as well as to
%(hopefully) make it easier to incorporate laser and camera functionality
%into other programs via the referenced functions.

% This program builds from code written by J.D. Rogers and Z. Poskin,
% as well as Automation2d.m by A. Radosevich.

%ZJS
%2018.7.1
%cleaned up slightly 
%2019.11.27

%%
%valid modes:
% fw
% bw
% abs

%test for valid mode:
if ~(strcmp(mode_flag,'fw')||strcmp(mode_flag,'bw')||strcmp(mode_flag,'abs'))
    disp('Error: not a valide mode.')
    return
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all

%assert for bw operation
% bw_flag=1;

%save the data out?
%savename = 'fixed_30um_onl_bbg'; % 0 will not save, set to a string to include in filename
%savename = 0; % 0 will not save, set to a string to include in filename

%galvo control:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('daqinitialized')
    s = daq.createSession('ni');
    addAnalogOutputChannel(s,'Dev1', 'ao0', 'Voltage');
    addAnalogOutputChannel(s,'Dev1', 'ao1', 'Voltage');
    daqinitialized = 1;
end
%center
scanparams.galvzero = [0 0];
s.outputSingleScan(scanparams.galvzero)

%laser set-up%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the handles for laser functions:
laser_session=laser_fns3;
%initialize, confusingly, we need the handle for the laser com port, called 'handles':
laser_handle=laser_session.laser_init();
pause(.5)
%set the default RF band to the visible band..
band=0;
%laser_session.laser_off(laser_handle)
%pause(.5)
laser_session.laser_RF_band(band,laser_handle)
pause(.5)
laser_session.laser_RF_band(band,laser_handle)%
pause(.5)
%turn laser on, i.e. call the laser_on fn in the session we set up, utilizing the handle:
laser_session.laser_RF_on(laser_handle)

%camera set-up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if bw_flag==1
if strcmp(mode_flag,'bw')
    %set up IDS or Quantalux camera..
    camera_selection='Quantalux';
    camera_handle = camera_handler(camera_selection);
    %camera_handle = camera_handler('Quantalux');
    %camera_handle = camera_handler('IDS');
else
    %set up ORCA camera:
    camera_selection='ORCA';
    camera_handle = camera_handler(camera_selection);
    %camera_handle = camera_handler('ORCA');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Galvo, Laser and Exposure time settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what colors do we want to look at?
%so all the settings are stored in a grid:
%color|RF level|exp time|laser power level

counter=1;
for colors=575:75:800
    run_settings(counter,1)=colors;
    counter=counter+1;
end

%set RF levels manually to adjust for different output
%if bw_flag==1
if strcmp(mode_flag,'bw')
    % backward: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %galvo voltages
    v_r=.7;
    voltages=[v_r v_r;v_r -v_r;-v_r -v_r;-v_r v_r];
    vc1=-.75:.25:-.5;
    vc2=.5:.25:.75
    voltages=[vc1' vc1';vc2' vc2'];
    
    %RF
    run_settings(1,2)=100;
    run_settings(2,2)=100;
    run_settings(3,2)=100;
    run_settings(4,2)=100;
    
    %exp time
    run_settings(:,3)=.05;
    %run_settings(:,3)=.025;
    run_settings(:,3)=.01;
    
    %laser power
    laser_power_lvl=100;
    laser_power_lvl=11;
    %laser_power_lvl=25;
    run_settings(:,4)=laser_power_lvl;
    %set laser power
    laser_session.laser_power(laser_power_lvl,laser_handle);
    
else
    % forward/abs: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RF
    %use same RF settings from abs to aid relating them... both fw and abs use
    %the same RF settings..
    
    %changed upon re-align for survey measurements
    run_settings(1,2)=9;
    run_settings(2,2)=9.7;
    run_settings(3,2)=4.5;
    run_settings(4,2)=8;
    
    %tuned 9-24-2019, is laser power changing too?
    run_settings(1,2)=7;
    run_settings(2,2)=12;
    run_settings(3,2)=8;
    run_settings(4,2)=12;
    
    %exp time, galvo positions -depends on fw or abs mode:
    if strcmp(mode_flag,'fw')
        %galvo voltages for fw
        start_voltage=.75;
        end_voltage=-.75;
        no_steps=4;
        step_voltage=(end_voltage-start_voltage)/no_steps;
        voltages_a=start_voltage:step_voltage:end_voltage;
        %to get toward 4pi coverage..
        voltages=[voltages_a' voltages_a';start_voltage -start_voltage;-start_voltage start_voltage];
        
        %exp time for fw
        run_settings(:,3)=.025;
        %run_settings(:,3)=.05;
        %run_settings(:,3)=.1;
        
    else
        %galvo positions for abs:
        voltages=[0 0];
        %exp time for abs:
        run_settings(:,3)=.002;
    end
    
    %laser power
    laser_power_lvl=11;
    run_settings(:,4)=laser_power_lvl;
    %set laser power
    laser_session.laser_power(laser_power_lvl,laser_handle)
end

[how_many columns]=size(run_settings);

%set exposure time: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_time=camera_handle.core.getExposure();%what are the units here?..
millisecs =1000*run_settings(1,3);
%max for TL camera is ~2000ms
camera_handle.core.setExposure(millisecs);% set exposure
exp_time=camera_handle.core.getExposure()%what are the units here?..
%image size
width=camera_handle.core.getImageWidth();
height=camera_handle.core.getImageHeight();

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Main Logic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hmm consider tightening up Laser commands for speed..
for color_n = 1:4%how_many
    %set laser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %logic to decide if we're changing RF band..set-up as visible..
    %if color_n==1, remains visible
    if (color_n==1 && run_settings((color_n),1)>=700)|| (color_n>1 && run_settings((color_n),1)>=700 && run_settings((color_n-1),1)<700)
        %change the band to IR, set color and RF
        disp('changing IR band')
        run_settings(color_n,:);
        band=1;
        laser_session.laser_RF_off(laser_handle);
        pause(.5)
        laser_session.laser_RF_band(band,laser_handle)
        pause(.5)
        laser_session.laser_RF_band(band,laser_handle)
        pause(.5)
        laser_settings=laser_session.laser_adj(run_settings(color_n,1),run_settings(color_n,2),laser_handle);
        pause(.5)
        laser_session.laser_RF_on(laser_handle)
        pause(.5)
        laser_session.laser_RF_on(laser_handle)
        pause(.5)
    else
        %set laser color and RF
        laser_settings=laser_session.laser_adj(run_settings(color_n,1),run_settings(color_n,2),laser_handle);
        pause(.5)
    end
    
    %decide on galvo positions
    counter=1; %for galvo positions:
    %{
    %if bw_flag==1
    if strcmp(mode_flag,'bw')
        %bw galvo positions
        v_r=.7;
        voltages=[v_r v_r;v_r -v_r;-v_r -v_r;-v_r v_r];
    elseif strcmp(mode_flag,'fw')
        start_voltage=.75;
        end_voltage=-.75;
        no_steps=4;
        step_voltage=(end_voltage-start_voltage)/no_steps;
        voltages_a=start_voltage:step_voltage:end_voltage;
        %to get toward 4pi coverage..
        voltages=[voltages_a' voltages_a';start_voltage -start_voltage;-start_voltage start_voltage];
        
    elseif strcmp(mode_flag,'abs')
    voltages=[0 0];
    end
    %}
    %{
    %set exposure time: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exp_time=camera_handle.core.getExposure();%what are the units here?..
    millisecs =1000*run_settings(1,3);
    %max for TL camera is ~2000ms
    camera_handle.core.setExposure(millisecs) ;% set exposure
    exp_time=camera_handle.core.getExposure()%what are the units here?..
    %image size
    width=camera_handle.core.getImageWidth();
    height=camera_handle.core.getImageHeight();
    %}
    %scan over positions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for angle_n=1:length(voltages)
    for angle_n=1:size(voltages,1)
        
        %move galvo, command via NI-6363
        s.outputSingleScan([voltages(angle_n,1) voltages(angle_n,2)])
        
        %average images?
        
        if strcmp(mode_flag,'abs')
            how_many_avgs=10;
        else
            how_many_avgs=1;
        end
        
        fail_count=0;
        
        for n=1:how_many_avgs
            disp([num2str(n) ' of ' num2str(how_many_avgs) ' to avg'])
            %{
                %background image
                for n=1:how_many_avgs
                    number=n
                    try  idsCAM_GRAB(hcam, I); %Stores the image in I.image
                        current_img=double(I.image');
                        img=img+current_img*1/how_many_avgs;
                    catch
                        fail_count=fail_count+1;
                    end
                end
                img=img*how_many_avgs/(how_many_avgs-fail_count);
                       
            %}
            % acquire %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fail_count=0;
            myflag = true;
            while myflag
                try
                    camera_handle.core.snapImage(); % take a shot
                    myflag=false;
                catch
                    fail_count=fail_count+1;
                    disp(strcat('snapImage() failed, consecutive fail count = ',num2str(fail_count)))
                    
                end
            end
            
            I1 = typecast(camera_handle.core.getImage(), 'uint16');
            %handles.frame.hide;
            I1 = double(I1);
            I1=reshape(I1,[width,height])';
            
            img=I1;
            
            %figure(15)
            %imagesc(img)
            %title(num2str(n))
            
            if strcmp(mode_flag,'abs')
                datacube(color_n,n,:,:)=img;
            end
            
            
        end %avgs at a postion
        
        %3D data cube %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~strcmp(mode_flag,'abs')
            datacube(color_n,angle_n,:,:)=img;
        end
        
        
        angle_n=angle_n+1;
        
        figure(15)
        imagesc(log(double(img)));
        hold on
        %circle(660,510,380)
        axis image
        if strcmp(mode_flag,'abs')
            title(['color ' num2str(run_settings(color_n,1)) ' max pixel val ' num2str(max(max(img)))])
            pause(1)
        else
            title(['color ' num2str(run_settings(color_n,1)) ' percent ' num2str(floor((angle_n-1)/length(voltages)*100))])
            pause(.5)
        end
        hold off
    end
    
    %beep;pause(1);beep;pause(1);beep   
    %beep off
    %Return to zero
    s.outputSingleScan(scanparams.galvzero)
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section to make a tag for keeping track of relevant info, settings, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag{1,1}='Sample Description'; %history:
tag{1,2}='0.52um silica spheres in glycerine and H2O for Bn(r) paper';
tag{1,2}='1.49um silica spheres in glycerine and H2O for Bn(r) paper, iris closed';
tag{1,2}='squirrel retina survey scan, ~80um pinhole..,bw w/ fw obj pulled away';
%tag{1,2}='1.54um polystyrene spheres in H2O, iris closed ~80um pinhole, trying to improve bw';
tag{1,2}='1um silica spheres in PDMS for Bn(r) paper, iris open';
tag{1,2}='liver and brain for Bn(r) paper, iris open, OD 3 for abs scat measurements';
%tag{1,2}='1.5 um ps in h2o, for camera comparison';
tag{1,2}='tendon scans to showcase asymmetry';
tag{1,2}='retina to try and nail down if averaging helps with us, ua';
tag{1,2}='lung to try and nail down if averaging helps with us, ua';
tag{1,2}='spheres, trying to nail down what is wrong with us, ua, hmm had to adjust laser RF';

tag{2,1}='Time of Experiment';
tag{2,2}=datetime;

%laser and exposure time settings:
tag{3,1}='Laser Wavelengths in nm/RF power percents/exposure times (ms)/laser powers (percents)';
tag{3,2}=run_settings;

%galvo settings:
tag{4,1}='galvo voltages';
tag{4,2}=voltages;

tag{5,1}='Polarization';
tag{5,2}='circular-ish';

%type of collection: fw, bw or abs
tag{6,1}='Forward, backward, or absolute:';
tag{6,2}=mode_flag;

%anything else
tag{7,1}='Camera:';
tag{7,2}=camera_selection;

%anything else
tag{8,1}='Additional Notes:';
tag{8,2}='nothing else to note';

%tag{7,2}='Interesting, seem to see hexagonal diffraction pattern from sphere packing?';
%tag{9,2}={'note: moving diagonally',...
%    'note2: black level, BLACK LEVEL OFFSET, gain and gamma correction all turned off in mex code',...
%    'note3: iris added before camera, oil in oil out, room lights off',...
%    'note4: monitor on, no wiggling, 10 avgs',...
%    'Xnote4: monitor off with manual wiggle,5 averages',...
%    'note5: bent fiber to get good unpolarized-looking output',...
%    'Xnote5: with electronic stage wiggling to smooth and kill speckle'};

%tag{6,2}={'Xnote: thinner sample: scotch tape spacer',...
%'note2: now using orca camera fw, TL bw, trying to optimize angles for actual 4pi coverage- more circular fw and bw, no on-center bw',...
%'Xnote3: wiht manual wiggle to eek out best contrast. '};


%%
%if savename, save(['quantalux/' datestr(now,'yyyymmddHHMMSS')  savename],'datacube','tag'); end

if savename, save(['gs_survey/' datestr(now,'yyyymmddHHMMSS')  savename mode_flag],'tag','datacube','-v7.3'); end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Disconnect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close galvo communication:
s.release; clear s;

%turn off RF:
laser_session.laser_RF_off(laser_handle)
pause(.5)
%close laser communication:
laser_session.laser_disconnect(laser_handle)

%close camera communication
camera_handle.core.unloadAllDevices(); % unload your device
%delete(instrfind)

end

