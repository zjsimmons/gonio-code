%functionalizing the laser stuff: using the functions in laser_fns.m:

%initialize the session for the laser, laser_fns returns handles to the
%various functions that do different laser operations, this hides all the
%bit command operations and should make the larger program easier to deal
%with.

%purpose of this is just a place to have/test the various laser fns. i.e.
%let's you do things like change the AOTF color and laser power from this
%script rather than with the NKT GUI. This is quite handy for applications
%where you might just want to spit out a laser color, etc.

%script also contains functionality to run the galvos. It's useful to have
%that by itself as well to facilitate debugging. 

%cleaned up a bit, 2019.11.27
%zjs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Laser Debug Functionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the handles for laser functions:
laser_session=laser_fns3

%initialize, confusingly, we need the handle for the laser com port, called
% 'handles':
handles=laser_session.laser_init();

%%
%turn it on, i.e. call the laser_on fn in the session we set up, utilizing the handle:
laser_session.laser_on(handles)
pause(.5)
laser_session.laser_RF_on(handles)

%%
%change the laser power
laser_settings=laser_session.laser_power(11,handles)

%%
%change the color, RF..
laser_settings=laser_session.laser_adj(575,5,handles);

%%
%loop through, rainbow light show:

pause(2)
for lambda=500:10:680
    lambda
    laser_settings=laser_session.laser_adj(lambda,10,handles);
    pause(1)
end

%{
lambda=530
for RF_p=1:.1:10
    RF_p;
    laser_settings=laser_session.laser_adj(lambda,RF_p,handles);
pause(1)
end
%}
%%
%note there are two output channels/two RF bands: vis and IR, you can only
%drive one at a time,

%change the RF band:
band=0;%this is vis channel, band=1 is IR

laser_session.laser_RF_off(handles)
pause(.5)
laser_session.laser_RF_band(band,handles)
%pause(1)
%laser_session.laser_RF_band(band,handles)


%hhmmm rf doesn't always seem to switch reliably, add pauses...? in
%in practice, I run this more than once when actually using this, still it
%seems like it will get hung up every once in a while.
%%
%turn RF on
laser_session.laser_RF_on(handles)

%%
%an IR color band selection
laser_settings=laser_session.laser_adj(725,100,handles)

%%
%turn laser off:
laser_session.laser_off(handles)

pause(.5)

%turn RF off:
laser_session.laser_RF_off(handles)

%%
%disconnect from laser
laser_session.laser_disconnect(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   other debug functonality: wiggle galvos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
if ~exist('daqinitialized')
    s = daq.createSession('ni');
    addAnalogOutputChannel(s,'Dev1', 'ao0', 'Voltage');
    addAnalogOutputChannel(s,'Dev1', 'ao1', 'Voltage');
    daqinitialized = 1;
end
%center
scanparams.galvzero = [0 0];
%scanparams.galvzero = [.45 .45];
s.outputSingleScan(scanparams.galvzero)

%%
%diagonal variations:
v_r=.7;
voltages=[v_r v_r;v_r -v_r;-v_r -v_r;-v_r v_r]

vc1=-.75:.25:-.5;
vc2=.5:.25:.75
voltages=[vc1' vc1';vc2' vc2'];


vc1=-.8:.1:.8;
%vc1=-.1:.001:.1;
voltages=[vc1' -vc1'];

for i = 1:1
    
    for angle_n=1:length(voltages)
        %move galvo, command via NI-6363
        %s.outputSingleScan([voltage -voltage])
        s.outputSingleScan([voltages(angle_n,1) voltages(angle_n,2)])
        voltages( angle_n,1)
        pause(.055)
        
    end
end
%}
scanparams.galvzero = [0 0];
%scanparams.galvzero = [-.8 -.8];

s.outputSingleScan(scanparams.galvzero)

%% circle scan
nv = 50;
maxv = 0.9899; % 0.8 is roughly the NA of the objective
maxv=.63
maxv=.8

for i=1:1
    for j=1:nv
        
        x = maxv*cos(j/nv*2*pi)% theta
        y = maxv*sin(j/nv*2*pi) % r gjm
        s.outputSingleScan([x y]);
        pause
        
    end
end

%scanparams.galvzero = [0 0];
%s.outputSingleScan(scanparams.galvzero)

%% do a scan quickly so we can watch
nv = 100;
maxv = 0.9; % 0.8 is roughly the NA of the objective
maxv = 0.1; % 0.8 is roughly the NA of the objective

for i=1:nv
    for j=1:nv
        s.outputSingleScan([i/nv-1/2 j/nv-1/2]*maxv);
        pause(.001)
        pause(.01)
    end
end


%% spiral scan
nv = 50;
maxv = 0.9; % 0.8 is roughly the NA of the objective
maxv=.1;
for i=1:nv
    for j=1:nv
        x = i/nv*maxv*cos(j/nv*2*pi);% theta
        y = i/nv*maxv*sin(j/nv*2*pi); % r
        s.outputSingleScan([x y]);
        pause(.001)
        pause(.02)
    end
end

scanparams.galvzero = [0 0];
s.outputSingleScan(scanparams.galvzero)


%%

%release galvos
s.release; clear s;
