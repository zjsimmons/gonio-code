%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Goniometry Survey Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function will translate the stage in steps and then run the 
%goniometry code for fw, bw, absolute detection at each step. This 
%functionality allows one to build up a linear tissue survey. 

%Note: This is absolute movement (in Zaber fn terms) because we don't want 
%to move it beyond the range of travel of our sample plate in the system.
%However: !IMPORTANT ISSUE: IF YOU UNPLUG THE ZABER STAGE, IT WILL GO HOME
%THE NEXT TIME YOU PLUG IT IN. THIS CAN RUN THE STAGE BEYOND THE CLEAR
%TRAVEL IN THE INSTRUMENT! 

%The main functionality that this script adds is translating the Zaber
%stage. The functions to do this are contained at the end of this program. 

%This code has 4 modes: 'fw_gonio', 'bw_gonio', 'transmission' and 'image'.
%The gonio and transmission modes all call the underlying gonio acquisition
%function 'gonio_fn_mm_ver3' (current ver) while 'image' mode opens the
%camera via this script and takes a sequence of images. 

%zjs
% cleaned up 2019.11.27

function gonio_survey_share(mode)

instrreset
base_filename='gs_';

%mode='transmission';
%%
%if in imaging mode, open the camera here, otherwise in gonio mode, the
%camera is opened by the gonio routine.

if strcmp(mode,'image')
    disp('opening the camera..')
    %  vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
    %  src = getselectedsource(vid);
    %  vid.FramesPerTrigger = 1;
    %  src.ExposureTime=0.05;
    %  src.ExposureTime=0.025;
    
    camera_selection='ORCA';
    camera_handle = camera_handler(camera_selection);
    
    millisecs =25;
    %max for TL camera is ~2000ms
    camera_handle.core.setExposure(millisecs) ;% set exposure
    exp_time=camera_handle.core.getExposure()%what are the units here?..
    %image size
    width=camera_handle.core.getImageWidth();
    height=camera_handle.core.getImageHeight();
    
end
%}

% % Open Zaber stage serial port
handles.zaber = serial('COM1','BaudRate',9600);
fopen(handles.zaber);

%disp('check in on zaber object:')
%handles.zaber

%this seems problematic to adjust speed, turned off -zjs
%{
%fwrite(handles.zaber,[1 42 160 15 0 0 ]) %set the speed to data of 4000, should be about max. Check manual, but should be like 3.7 mm/s
%hmmm this seems to be causing problems, seems to be getting stuck..
%slower...1000 instead of 4000
%fwrite(handles.zaber,[1 42 232 3 0 0 ]) %set the speed to data of 4000, should be about max. Check manual, but should be like 3.7 mm/s
%hmm this doesn't seem to move as expected..work over the zaber fn parts..
%}

%let's set the step size too
step_size=32;
fwrite(handles.zaber,[0 37 step_size 0 0 0 ])%that means about 32 steps per microstep..

%how much total to translate the stage:

%after rebuild, larger travel.. may be better to set this via writing max and min travel
%values, but this didn't work through the zaber console and I don't want to
%go to the trouble trying to implement it here. As a result:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       WARNING
%       Do not make the offset less than 7500 or the travel more than 6mm
%       as this would allow for crashing the little sample holder into the
%       uprights in the goniometer set-up.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%used with lung:
max_travel_mm=1;
offset=8500;

%larger range- this is totally fine
max_travel_mm=4;
max_travel_mm=1;

offset=8500;
offset=7500;

%go home
%go_home(handles);
%fwrite(handles.zaber,[0 1 0 0 0 0 ])
%pause(10)

%move to starting position
move_stage_abs(offset,handles)
pause(5)

no_positions=25;%no of locations to sample
no_positions=30;
no_positions=15;
%no_positions=10;

%no_positions=1;max_travel_mm=1.2;max_travel_mm=.1;
%no_positions=3;max_travel_mm=0;max_travel_mm=.1;

phase_offset=-33;%microns as well
phase_offset=0;

offset=offset+phase_offset;
for n=1:no_positions
    
    %I1 = take_im_and_update(handles);
    %absolute distance to move in microns
    distance=offset+(n)/no_positions*1000*max_travel_mm
    %position in mm
    pos(n+1)=distance/1000;
    move_stage_abs(distance,handles)
    pause(.5)
    
    %this would be where you call the gonio acq script, for each position, need
    %a version that takes the fw/bw option as an input..
    
    if strcmp(mode,'fw_gonio')
        disp('fw goniometry mode...')
        
        %
        disp('forward measurement proceeding')
        %call for fw
        %gonio_fn(strcat('1pt49_um_fw',num2str(n)),0);
        %gonio_fn(strcat('sqr_retina_seq_fw',num2str(n)),0);
        gonio_fn_mm_ver3(strcat(base_filename,num2str(n)),'fw');
        
        pause(1)%rationale that this is waiting a beat to let the stage settle
        
        %wait for button press- do we take up the fw objective?
        % disp('ready for BW?, press a key')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
        % pause;     
        %}
     
    elseif strcmp(mode,'bw_gonio')
        disp('backward measurement proceeding')
        disp('XXXXX dont forget to take out the ND filters XXXX')
        
        %call for bw
        %gonio_fn(strcat('1pt49_um_bw',num2str(n)),1);
        %  gonio_fn(strcat('sqr_retina_seq_bw',num2str(n)),1);
        gonio_fn_mm_ver3(strcat(base_filename,num2str(n)),'bw');
        
        %}     
        pause(1)
        %wait for button press- do we take up the fw objective?
        %disp('ready for FW again?, press a key')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
        %pause;
        %disp('forward measurement proceeding')
          
    elseif strcmp(mode,'transmission')
        disp('absolute transmission measurement')
        disp('Did you remember to put in the ND filter?')
        %add logic for calling the abs measurement routine
        %gonio_abs_fn(strcat('sqr_retina_abs',num2str(n)));
        %gonio_fn_mm_ver2(strcat('sqr_retina_',num2str(n)),'abs');
        gonio_fn_mm_ver3(strcat(base_filename,num2str(n)),'abs');
        %    gonio_fn_mm_ver3(strcat('black_',num2str(n)),'abs');
        
        pause(1)
       
        handles.zaber
        
    elseif strcmp(mode,'image')
        disp('imaging mode...')
               
        try
            camera_handle.core.snapImage(); % take a shot
        catch
            disp('problem with image acquisition, snapImage() failed')
            pause
            %do something here in case it fails? need to test.. naively try
            %again..
            camera_handle.core.snapImage(); % take a shot            
        end
        
        I1 = typecast(camera_handle.core.getImage(), 'uint16');
        %handles.frame.hide;
        I1 = double(I1);
        I1=reshape(I1,[width,height])';     
        img=I1;
        img_seq{n+1}=img;
                
        figure(15)
        imagesc(img)
        hold on
        title(strcat('squirrel retina seq ',num2str(n)))     
        %can display (manually adjust) approx illuminated location with
        %iris closed. (if desired)
        circx=1300;circy=800;circr=400;
        %circle(circx,circy,circr)%w/ this set-up, 0.1657 um/pix, 360pix=60um
        hold off
        
        pause(1)
    else
        disp('problem with mode selction')
    end
    
    %  disp('check in on zaber object:')
    %  handles.zaber
    
end% scanning over positions

%go back
move_stage_abs(offset,handles)
%close stage
delete(handles.zaber)

%if in imaging mode, save the sequence..
if strcmp(mode,'image')
    savename=strcat(base_filename,'image_seq_open_iris_13lgs');
    % savename=0;
    
    % Section to make a tag for keeping track of relevant info, settings, etc.
    tag{1,1}='Sample Description';
    tag{1,2}='scan imaging sequence- 13LGS retina';
    
    tag{2,1}='Time of Experiment';
    tag{2,2}=datetime;
    
    tag{3,1}='Illumination';
    tag{3,2}='flashlight,  100um spot';
    
    tag{4,1}='approx illumination location';
    tag{4,2}=[circx circy circr];
    
    %polarization
    %tag{5,1}='Positions';
    %tag{4,2}=pos;
    
    %anything else
    tag{5,1}='Additional Notes:';
    tag{5,1}='this set-up, 0.1657 um/pix, 360pix=60um';
    
    %
    if savename, save(['gs_survey/' datestr(now,'yyyymmddHHMMSS')  savename],'img_seq','tag'); end
    
end

end

%%
%move the stage home
function []=go_home(handles)
fwrite(handles.zaber,[0 1 0 0 0 0 ])
end

%%
%move the stage to an absolute position
function []=move_stage_abs(distance,handles)
step_size=32;

%desired_distance=500;%microns
desired_distance=distance;
steps=round(desired_distance/.09921875*step_size/128);%note: apparently for this device (LS28E), the default step is 128 when converting to abs position.
%convert to bytes:
cmd_data=steps;

%'Handles negative data- incorporating from zaber
if cmd_data < 0
    cmd_data = 256^4 + cmd_data;
end
%the four data bytes:
cmd_byte_6 = floor(cmd_data/256^3);
cmd_data = cmd_data-256^3*cmd_byte_6;
cmd_byte_5 = floor(cmd_data/256^2);
cmd_data = cmd_data-256^2*cmd_byte_5;
cmd_byte_4 = floor(cmd_data/256);
cmd_data = cmd_data-256*cmd_byte_4;
cmd_byte_3=cmd_data;

fwrite(handles.zaber,[0 20 cmd_byte_3 cmd_byte_4 cmd_byte_5 cmd_byte_6])
%disp('stage moved')
end


