function [handles] = camera_handler(camera_type)


%camera_communication_sample.m
%
%So this is a bit of code that contains how to communicate with a camera
%via micromanager. Goal is for this to be about the minimal piece of code
%needed to aid in integrating imaging via micromanager into other
%programs. Mostly excised from Radosevich Automation2d.m, but also modified 
%to get it to work with MM 2.0 and to incorporate different cameras (IDS 
%and ORCA flash, and most recently Quantalux).

%note: this will require micromanager and adapters for specific cameras. If
%you can get the camera working in the micromanager gui, should work via
%this method as well.  

%zjs 2018.6.28



%%
%set-up

addpath(genpath('C:\Program Files\Micro-Manager-2.0beta'))  

import mmcorej.* % load the micro-manager library
%import MMCOREJ.* % load the micro-manager library

handles.core = CMMCore; %instantiate a camera object
handles.core.reset;
handles.core.unloadAllDevices(); % unload your device
%handles.core.loadDevice('Camera1', 'PrincetonInstruments', 'Camera-1'); % load your camera.

%choose camera:
%camera_type='IDS';
%camera_type='Quantalux';
%camera_type='ORCA';
%expand to other cameras??

if strcmp(camera_type,'ORCA')
    handles.core.loadDevice('Camera1', 'HamamatsuHam', 'HamamatsuHam_DCAM'); % load your camera.
elseif strcmp(camera_type,'Quantalux')
    handles.core.loadDevice('Camera1', 'TSI', 'TSICam');%load your camera.
elseif strcmp(camera_type,'IDS')
%original recipe:
   % handles.core.loadDevice('Camera1', 'IDS_uEye_original', 'IDS uEye'); % load your camera.
%modified version that incorporates changes..
    %handles.core.loadDevice('Camera1', 'IDS_uEye', 'IDS uEye ver2hp'); % load your camera.
    handles.core.loadDevice('Camera1', 'IDS_uEye_modified', 'IDS uEye ver2hp'); % load your camera.
else
    disp('problem with camera selction')
end

% handles.core.loadDevice('Camera1', 'PVCAM', 'Camera-1'); % load your camera.
%See if camera loads, otherwise throw an error
try
    %so this doesn't seem to work, it can execute, but subsequent calls
    %won't work..zjs
    %handles.core.initializeDevice('Camera1'); %initialize my device
    handles.core.initializeAllDevices();%hmm the more general call does work
catch
    errordlg('Camera Didn''t Initiate''')
    return;
end

%don't need to set the temp for the Orca..zjs
    %handles.core.setProperty('Camera1','CCDTemperatureSetPoint','-75');
    %set(handles.actual_temp,'String',char(handles.core.getProperty('Camera1','CCDTemperature')))

% % % % check property names and values for camera
 k = handles.core.getDevicePropertyNames('Camera1');
 for i=0:k.size-1;
 prop(i+1) = k.get(i);
 val(i+1) = handles.core.getProperty('Camera1', prop(i+1));
 disp([char(prop(i+1)) ':  ' char(val(i+1))])
 %pause(.4)
 end

 
disp('now change settings:')
%hmm looks like we can set the frame rate thus, and all the other available
%properties as well: 
if strcmp(camera_type,'IDS')
%set frame rate:
propFrameRate = handles.core.getProperty('Camera1', 'Frame Rate')
handles.core.setProperty('Camera1', 'Frame Rate', '.5');
propFrameRate = handles.core.getProperty('Camera1', 'Frame Rate')

%set gain, usually off: 
propGain = handles.core.getProperty('Camera1', 'Gain')
%handles.core.setProperty('Camera1', 'Gain', '50');
propGain = handles.core.getProperty('Camera1', 'Gain')

%OK, let's set our new properties that we have access to: 
%note: don't have access to all of these properties w/ original adapter)
%set offset (black level), 255 for now: 
propBlackLevel = handles.core.getProperty('Camera1', 'Offset')
handles.core.setProperty('Camera1', 'Offset', '255');
propBlackLevel = handles.core.getProperty('Camera1', 'Offset')

%set shutter mode: 
propShutterMode = handles.core.getProperty('Camera1', 'Sensor Shutter Mode')
handles.core.setProperty('Camera1', 'Sensor Shutter Mode', '1- Rolling shutter');
%handles.core.setProperty('Camera1', 'Sensor Shutter Mode', '2- Global shutter');
%handles.core.setProperty('Camera1', 'Sensor Shutter Mode', '1- Rolling shutter');
propShutterMode = handles.core.getProperty('Camera1', 'Sensor Shutter Mode')

%set log mode
propLogMode = handles.core.getProperty('Camera1', 'Sensor Log Mode')
handles.core.setProperty('Camera1', 'Sensor Log Mode', '1- log mode off');
propLogMode = handles.core.getProperty('Camera1', 'Sensor Log Mode')

%set hot pixel correction
propHPCMode = handles.core.getProperty('Camera1', 'Hotpixel correction (sensor)')
handles.core.setProperty('Camera1', 'Hotpixel correction (sensor)', 'sensor hot pixel correction on');
%handles.core.setProperty('Camera1', 'Hotpixel correction (sensor)', 'sensor hot pixel correction off');
propHPCMode = handles.core.getProperty('Camera1', 'Hotpixel correction (sensor)')

%hmm pixel clock? 
propPixelClock = handles.core.getProperty('Camera1', 'Pixel Clock')
handles.core.setProperty('Camera1', 'Pixel Clock', '7');
%handles.core.setProperty('Camera1', 'Hotpixel correction (sensor)', 'sensor hot pixel correction off');
propPixelClock = handles.core.getProperty('Camera1', 'Pixel Clock')

elseif strcmp(camera_type,'ORCA')
%set conversion factor coefficient 
propCFactor = handles.core.getProperty('Camera1', 'CONVERSION FACTOR COEFF')
handles.core.setProperty('Camera1', 'CONVERSION FACTOR COEFF', '1');
propCFactor = handles.core.getProperty('Camera1', 'CONVERSION FACTOR COEFF')

%set conversion factor offset
propCFOffset = handles.core.getProperty('Camera1', 'CONVERSION FACTOR OFFSET')
handles.core.setProperty('Camera1', 'CONVERSION FACTOR OFFSET', '0');
propCFOffset = handles.core.getProperty('Camera1', 'CONVERSION FACTOR OFFSET')

%hmm appears you can't write these vals? interesting, looks like these are
%fixed and writing them doesn't change anything.. 

elseif strcmp(camera_type,'Quantalux')
    disp('quantalux hello')
%set hotpixel correction.. 
propHPmode = handles.core.getProperty('Camera1', 'HotPixel')
%handles.core.setProperty('Camera1', 'HotPixel', 'On');
propHPmode = handles.core.getProperty('Camera1', 'HotPixel')

end

end

%re: Quantalux, really not a lot of parameters to set, you have in TSI
%anyway, hot pixel correction and EEP, and of course exposure time. 

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%everything above is camera specific, below is not:


%%
%Take an image:

width=handles.core.getImageWidth();
height=handles.core.getImageHeight();
    
%set exposure time: 
%so here's functions for getting and setting the exposure time: 
exp_time=handles.core.getExposure()%what are the units here?..
millisecs =10 
%max for TL camera is ~2000ms
handles.core.setExposure(millisecs) % set exposure
exp_time=handles.core.getExposure()%what are the units here?..

%take image:
handles.core.snapImage(); % take a shot
I1 = typecast(handles.core.getImage(), 'uint16');
%handles.frame.hide;
I1 = double(I1);
I1=reshape(I1,[width,height])';
    
figure(12)
imagesc(I1)
axis image;axis off;

%%
%take images:
        max_reps=10
        reps=1;
        while reps <= max_reps
            handles.core.snapImage(); % take a shot
            temp = typecast(handles.core.getImage(), 'uint16');
            temp = double(temp);
            temp=reshape(temp,[width,height])';
            figure(12)
            imagesc(temp)
            title(num2str(reps))
            axis image;axis off;           
            % I1 = I1+temp;
            reps = reps+1;
        end


%%
%close the camera
handles.core.unloadAllDevices(); % unload your device
delete(instrfind)

%}