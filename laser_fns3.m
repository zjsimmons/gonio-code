function fh = laser_fns3()

%This is a functionalized laser control program, built from  
%LaserControl.m and setlaserwavelength.m by JDR and ZCP (which were adapted 
%from Andrew Radosevich's automation2d.m script)

%This program contains all the conversions to telegram messages for different 
%operations as functions. By making the different operations functions 
%whose handles are returned by this, should make it much cleaner to change 
%the laser settings via calling the hadles in other programs.. and so to 
%integrate it into the data acquisition progogram, etc. 

% note re: addresses: '015' is the laser, '016' is RF driver, '017' is
% superK select, 

%zjs, 2017.3.27, updated 2018.1.5

%updated 2018.7.12, 'special character' bug fix, see: tele_check.m

fh.laser_init=@initialize_laser_com;
fh.laser_disconnect=@close_laser_com;
fh.laser_on=@turn_on_laser;
fh.laser_RF_on=@turn_on_RF;
fh.laser_off=@turn_off_laser;
fh.laser_adj=@set_laser;
fh.laser_RF_band=@set_RF_band;
fh.laser_RF_off=@turn_off_RF;
fh.laser_power=@set_laser_power;


end

%% fn to initialize NKT laser module.
function [handles]=initialize_laser_com()
handles.laser = serial('COM3');
set(handles.laser,'BaudRate',115200,'StopBits',1,'Terminator','');
set(handles.laser,'DataTerminalReady','off');%must do this first for some reason in order to read.
fopen(handles.laser);
set(handles.laser,'DataTerminalReady','on');%must do this first for some reason in order to read.
end

%% fn close communication with the NKT laser
function close_laser_com(handles)
fclose(handles.laser);
delete(handles.laser);
clear handles;
end

%%
%{
% check power
% '013' is start telegram. '015' is SuperK destination module. '166' is
% computer address. '004' is read ('005' is write). '055' is register for
% SuperK power.
string2 = ['013'; '015';'166';'004';'055'];
OldCRC = 0;
for i = 2:length(string2);
    OldCRC = crc16_add1(str2num(string2(i,:)),OldCRC);
end
temp = dec2hex(OldCRC);
temp = [num2str(zeros(1,4-length(temp))) temp];
seven = '000';
seven(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
eight = '000';
eight(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
string2 = [string2 ;seven ;eight;'010'];

%read out the current laser power:
fwrite(handles.laser, uint8(str2num(string2)),'uint8');
pause(0.15);
back = fread(handles.laser,get(handles.laser,'BytesAvailable'));
power = back(6,:) + back(7,:)*256
%}
%% fn Turn on emission 
function[]=turn_on_laser(handles)
% '048' is SuperK emission register. '003' is data byte for emission on.
string = ['013';'015';'166';'005';'048';'003';'118';'016';'010'];%Command to turn on emission
fwrite(handles.laser, uint8(str2num(string)),'uint8');
end

%% fn Turn on RF power
function[]=turn_on_RF(handles)
% '048' is also RF driver power register. '001' is data byte for RF power on
string = ['013';'016';'166';'005';'048';'001';'055';'241';'010'];%Command to turn on RF power
fwrite(handles.laser, uint8(str2num(string)),'uint8');
end

%% change RF band/crystal to IR:
% '052' is RF band selection register. '001' is data byte for turning on 2nd RF band on.
% note: this is sent to the superK select module, address: '017'. 
%string = ['013';'015';'166';'005';'052';'001';'118';'016';'010'];%Command to turn on emission
%fwrite(handles.laser, uint8(str2num(string)),'uint8');
%
function[]=set_RF_band(band,handles)

IR_on=band;

if IR_on==1
    switch_string='001';
elseif IR_on==0
    switch_string='000';
else
    disp('error, problem with RF band selection');
end

string2 = ['013'; '017';'166';'005';'052';switch_string];
OldCRC = 0;
for i = 2:length(string2)
    OldCRC = crc16_add1(str2num(string2(i,:)),OldCRC);
end
temp = dec2hex(OldCRC);
temp = [num2str(zeros(1,4-length(temp))) temp];
seven = '000';
seven(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
eight = '000';
eight(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
string2 = [string2 ;seven ;eight;'010'];
fwrite(handles.laser, uint8(str2num(string2)),'uint8');
pause(.15)
%this just spits out what's coming back from the device..
while get(handles.laser,'BytesAvailable')>0; buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'));
end

end
%}
%% %fn to set the laser color and power
function[laser_settings]=set_laser(color,RF_level,handles,varargin)

%first time through, initialize all noted laser settings to zero:
if ~exist('laser_settings')
    laser_settings.color=0;
    laser_settings.RF_level=0;
    laser_settings.RF_band=0;
    laser_settings.power_level=0;
end

%initially tried to incorporate automatic switching of RF bands; this
%didn't work, perhaps a matter of having adequate pause between operations?
%In any case, only need to do this infrequently, so will be done more
%manually with an RF band selection function. 
%{
%determine what the band is... if switching, turn off and switch band..
if laser_settings.color<700 && color>700
    %switch to IR
    disp('switch to IR')
   
    %turn off the RF
    % '048' is also RF driver power register. '000' is data byte for RF power off
    string = ['013';'016';'166';'005';'048';'000';'039';'208';'010'];%Command to turn off RF power
    fwrite(handles.laser, uint8(str2num(string)),'uint8');
    pause(0.15);%ok, in order to get both commands through, color and amplitude, seems to need a pause, -zjs

    
    %switch the RF
    switch_string='001';
    string2 = ['013'; '017';'166';'005';'052';switch_string];
    OldCRC = 0;
    for i = 2:length(string2);
        OldCRC = crc16_add1(str2num(string2(i,:)),OldCRC);
    end
    temp = dec2hex(OldCRC);
    temp = [num2str(zeros(1,4-length(temp))) temp];
    seven = '000';
    seven(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
    eight = '000';
    eight(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
    string2 = [string2 ;seven ;eight;'010'];
    %write
    fwrite(handles.laser, uint8(str2num(string2)),'uint8');
    %this just spits out what's coming back from the device..
    while get(handles.laser,'BytesAvailable')>0;buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'))
        end;
    pause(0.15);
    %have to write it twice for some reason
    fwrite(handles.laser, uint8(str2num(string2)),'uint8');
    %this just spits out what's coming back from the device..
    while get(handles.laser,'BytesAvailable')>0;buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'))
        end;
    pause(0.15);
    %turn the RF back on:     
    % '048' is also RF driver power register. '001' is data byte for RF power on
    string = ['013';'016';'166';'005';'048';'001';'055';'241';'010'];%Command to turn on RF power
    fwrite(handles.laser, uint8(str2num(string)),'uint8');
    pause(0.15);%ok, in order to get both commands through, color and amplitude, seems to need a pause, -zjs

    laser_settings.RF_band=1;
    
elseif laser_settings.color>700 && color<700
    disp('switch to vis')  
    
    %turn off the RF
    % '048' is also RF driver power register. '000' is data byte for RF power off
    string = ['013';'016';'166';'005';'048';'000';'039';'208';'010'];%Command to turn off RF power
    fwrite(handles.laser, uint8(str2num(string)),'uint8');
    pause(0.15);%ok, in order to get both commands through, color and amplitude, seems to need a pause, -zjs

    
    %switch the RF
    switch_string='000';
    string2 = ['013'; '017';'166';'005';'052';switch_string];
    OldCRC = 0;
    for i = 2:length(string2);
        OldCRC = crc16_add1(str2num(string2(i,:)),OldCRC);
    end
    temp = dec2hex(OldCRC);
    temp = [num2str(zeros(1,4-length(temp))) temp];
    seven = '000';
    seven(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
    eight = '000';
    eight(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
    string2 = [string2 ;seven ;eight;'010'];
    %write
    fwrite(handles.laser, uint8(str2num(string2)),'uint8');
    %this just spits out what's coming back from the device..
    while get(handles.laser,'BytesAvailable')>0;buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'))
        end;
    pause(0.15);
    %have to write it twice for some reason
    fwrite(handles.laser, uint8(str2num(string2)),'uint8');
    %this just spits out what's coming back from the device..
    while get(handles.laser,'BytesAvailable')>0;buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'))
        end;
    pause(0.15);
    %turn the RF back on:     
    % '048' is also RF driver power register. '001' is data byte for RF power on
    string = ['013';'016';'166';'005';'048';'001';'055';'241';'010'];%Command to turn on RF power
    fwrite(handles.laser, uint8(str2num(string)),'uint8');
    pause(0.15);%ok, in order to get both commands through, color and amplitude, seems to need a pause, -zjs

    
    laser_settings.RF_band=0;
else
    disp('no change in RF band')  
end %RF band switch..
%}

if laser_settings.color~=color
% Determine string for channel selected
% Wavelength channels run from register 144-151 for channels 1-8
channel=1;%we're only doing one at a time..
wavechannel = num2str(channel+143);

%%%%%%%%%%%%%%%% Create telegram for setting wavelength%%%%%%%%%%%%%%%%%%%%
% '013' is start message. '016' is RF driver address. '166' is source
% address (computer). '005' is write ('004 is read).
string = ['013';'016';'166';'005'];

Wavelength = color*1000; % convert wavelength from nm to pm to read in
% Convert wavelength from dec to hex, then rearrange bytes from LSB to MSB
% to read into the telegram
Wavelength = dec2hex(Wavelength,8);
six = '000'; temp = num2str(hex2dec(num2str(Wavelength(7:8)))); six((3-length(temp)+1):end) = temp;
seven = '000'; temp = num2str(hex2dec(num2str(Wavelength(5:6)))); seven((3-length(temp)+1):end) = temp;
eight = '000'; temp = num2str(hex2dec(num2str(Wavelength(3:4)))); eight((3-length(temp)+1):end) = temp;
nine = '000'; temp = num2str(hex2dec(num2str(Wavelength(1:2)))); nine((3-length(temp)+1):end) = temp;

% Calculate CRC values with updated telegram string
string = [string;wavechannel;six;seven;eight;nine];
OldCRC = 0;
for i = 2:length(string)
    OldCRC = crc16_add1(str2num(string(i,:)),OldCRC);
end
temp = dec2hex(OldCRC);
temp = [num2str(zeros(1,4-length(temp))) temp];
ten = '000';
ten(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
eleven = '000';
eleven(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
stringwave = [string;ten;eleven;'010'];
%check string for 'special characters': 
stringwave=tele_check(stringwave);
%write: color
fwrite(handles.laser, uint8(str2num(stringwave)),'uint8');
pause(0.15);%ok, in order to get both commands through, color and amplitude, seems to need a pause, -zjs

laser_settings.color=color;
end

if laser_settings.RF_level~=RF_level
%%%%%%%%%%%%%%%% Create telegram for setting amplitude%%%%%%%%%%%%%%%%%%%%
string = ['013';'016';'166';'005'];
% Wavelength amplitude channels run from register 176-183 for channels 1-8
five = num2str(channel+175);
Amplitude = RF_level*10; %Amplitude is read in as tenths of a percent
Amplitude = cast(Amplitude,'uint16');
% Convert amplitude from dec to hex, then rearrange bytes from LSB to MSB
% to read into the telegram
Amplitude = dec2hex(Amplitude,4);
six = '000'; temp = num2str(hex2dec(num2str(Amplitude(3:4)))); six((3-length(temp)+1):end) = temp;
seven = '000'; temp = num2str(hex2dec(num2str(Amplitude(1:2)))); seven((3-length(temp)+1):end) = temp;
% Calculate CRC values with updated telegram string
string = [string;five;six;seven];
OldCRC = 0;
for i = 2:length(string)
    OldCRC = crc16_add1(str2num(string(i,:)),OldCRC);
end
temp = dec2hex(OldCRC);
temp = [num2str(zeros(1,4-length(temp))) temp];
eight = '000';
eight(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
nine = '000';
nine(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
stringamp = [string;eight;nine;'010'];
%check for 'special characters'
stringamp=tele_check(stringamp);

%write amplitude:
fwrite(handles.laser, uint8(str2num(stringamp)),'uint8');
pause(.15)

laser_settings.RF_level=RF_level;
end



%setlaserwavelength(1,600,100,handles);
%setlaserwavelength(1,color,rf_power,handles);
%pause(.5)
%while get(handles.laser,'BytesAvailable')>0;buffer_contents=fread(handles.laser,get(handles.laser,'BytesAvailable'))
%end;

end
%%
%IR?
%setlaserwavelength(1,775,100,handles);
%pause(.5)
%% Set wavelength and amplitude (AOTF power)
%{
channel = 1;
amplitude = 1;

a=tic
for j =1:3:10
    for i = 1:10:100
       color=549+i
       power=10*j
    setlaserwavelength(channel,color,power,handles)
    pause(.5)
    end
end

toc(a)
%}

%% Adjust laser power.. 
function[laser_settings]=set_laser_power(power_level,handles,varargin)

%apparently register is 37h.. 

%if laser_settings.power_level~=power_level
    
if ~exist('laser_settings')
    laser_settings.color=0;
    laser_settings.RF_level=0;
    laser_settings.RF_band=0;
    laser_settings.power_level=0;
end    

if laser_settings.power_level~=power_level

%%%%%%%%%%%%%%%% Create telegram for setting amplitude%%%%%%%%%%%%%%%%%%%%
% '013' is start message. '015' is superK '166' is source
% address (computer). '005' is write ('004 is read).
string = ['013';'015';'166';'005'];
% Wavelength amplitude channels run from register 176-183 for channels 1-8
five = '055';
Amplitude = power_level*10; %Amplitude is read in as tenths of a percent
% Convert amplitude from dec to hex, then rearrange bytes from LSB to MSB
% to read into the telegram
Amplitude = dec2hex(Amplitude,4);
six = '000'; temp = num2str(hex2dec(num2str(Amplitude(3:4)))); six((3-length(temp)+1):end) = temp;
seven = '000'; temp = num2str(hex2dec(num2str(Amplitude(1:2)))); seven((3-length(temp)+1):end) = temp;

% Calculate CRC values with updated telegram string
string = [string;five;six;seven];
OldCRC = 0;
for i = 2:length(string)
    OldCRC = crc16_add1(str2num(string(i,:)),OldCRC);
end
temp = dec2hex(OldCRC);
temp = [num2str(zeros(1,4-length(temp))) temp];
eight = '000';
eight(3-length(num2str(hex2dec(temp(1:2))))+1:end) = num2str(hex2dec(temp(1:2)));
nine = '000';
nine(3-length(num2str(hex2dec(temp(3:4))))+1:end) = num2str(hex2dec(temp(3:4)));
stringamp = [string;eight;nine;'010'];
%check for 'special characters'
stringamp=tele_check(stringamp);

%write amplitude:
fwrite(handles.laser, uint8(str2num(stringamp)),'uint8');
pause(1)

laser_settings.power_level=power_level;
end


end




%% Turn off emission
function[]=turn_off_laser(handles)
% '048' is SuperK emission register. '000' is data byte for emission off.
string = ['013';'015';'166';'005';'048';'000';'070';'115';'010'];%Command to turn off emission
fwrite(handles.laser, uint8(str2num(string)),'uint8');
end

%% Turn off RF power

function[]=turn_off_RF(handles)

% '048' is also RF driver power register. '000' is data byte for RF power off
string = ['013';'016';'166';'005';'048';'000';'039';'208';'010'];%Command to turn off RF power
fwrite(handles.laser, uint8(str2num(string)),'uint8');

pause(.5)

end

