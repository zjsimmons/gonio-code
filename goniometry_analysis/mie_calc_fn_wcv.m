function [norm_factor,varargout]=mie_calc_fn_wcv(lambda,diameter,cv,sphere_n,bulk_n,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     MIE calculation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DESCRIPTION: This function calculates Mie scattering curves as a function 
%theta. This is for comparison to experimentally-derived Mie data.
%Currently it outputs the curve for 3 polarizations: | _ and equal, i.e.
%both. There are a variety of output plot options, see code. A log plot
%seems to be one of the more useful ones. 

%note: in terms of the math, this program does the actual calcuation, 
%finding the Bessel fn expansion solution up to a given no of terms, set at
%100 by default. The two called fns: a_n_b_n_fn and pi_tau_fn find the
%Bessel expansion coefficients and angle-dependent terms respectively. See
%Bohren and Huffman ch. 4 for more info, specifically 4.46, 4.56, 4.74,etc.  

%note: this contains a variety of estimation/extension logic but this 
%current function version does not utilize it. The logic is in the 
%commented section at the bottom (majority of file) and may be incorporated
%in a ver2 of the function implementation. 

%INPUTS: wavelength of interest (lambda) in um, and sphere diameter in um 

%refractive indicies are currently hard-coded but could be made into inputs
%as well. 

%OUTPUTS: Displays Mie scattering curve. 
%norm_factor - for normalization purposes when working with experimental
%traces.. 
%varargout{1}=thetas*180/pi;
%varargout{2}=intensity with theta angle;



%FUNCTIONS CALLED: 
%a_n_b_n_fn(m,x,n)
%pi_tau_fn

%zjs, written fall 2015, cleaned up 10-24-2016

%major revisions 2018: added material options, cv input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pol_to_plot=varargin{1}; %output display option

%%
%some parameters: 
d=diameter;
a=d/2;

%for reference:
%N_pdms=1.412;%from zach 1
%N_silica=1.4589;%from zach 1
%N_H2O=1.3338;
%N_polystyrene=1.5931;
%N_glycerine=1.4729;
%N_H2O=1.337;

% Better, spectral:
% H2O
C1 = 5.684027565E-1; C2 = 5.101829712E-3; C3 = 1.726177391E-1; 
C4 = 1.821153936E-2; C5 = 2.086189578E-2; C6 = 2.620722293E-2; 
C7 = 1.130748688E-1; C8 = 1.069792721E1;
N_H2O = sqrt(1 + C1*lambda.^2./(lambda.^2-C2) +...
    C3*lambda.^2./(lambda.^2-C4) + C5*lambda.^2/(lambda.^2-C6) +...
    C7*lambda.^2/(lambda.^2-C8));% clear C1 C2 C3 C4 C5 C6 C7 C8

% glycerine
N_glycerine = 1.45797 + .00598*lambda.^-2 -.00036 *lambda^-4; 

% polystyrene
N_polystyrene = 1.5663 + 7.85e-3*lambda.^-2 + 3.34e-4*lambda.^-4;

%PDMS Sylgard 184 see F. Schneider 2009
%npdms=1.4;
B1=1.0093; 
C1=13185; %nm^2
C1=C1/1000^2; %um^2
N_pdms=sqrt(1+B1*lambda.^2/(lambda.^2-C1));

% silica
ns1=.6961663*lambda^2/(lambda^2-.0684043^2);
ns2=.4079426*lambda^2/(lambda^2-.1162414^2);
ns3=.8974794*lambda^2/(lambda^2-9.896161^2);
N_silica=sqrt(ns1+ns2+ns3+1);

%sphere index options: 'silica', 'polystyrene'
if strcmp(sphere_n,'silica')
n_1=N_silica;
elseif strcmp(sphere_n,'polystyrene')
n_1=N_polystyrene;
else
disp('Problem: invalid sphere index selection')
end

%sphere index options: 'H2O', 'PDMS', 'glycerine'
if strcmp(bulk_n,'H2O')
n_0=N_H2O;
elseif strcmp(bulk_n,'PDMS')
n_0=N_pdms;
elseif strcmp(bulk_n,'glycerine')
n_0=N_glycerine;
else
disp('Problem: invalid bulk index selection')
end


n_0=n_0*1;

%%
%sphere size distribution gives weighting.. 

if cv~=0

sigma=cv*diameter;
x_min=diameter-3*sigma;
x_max=diameter+3*sigma;

no_samples=50;%keep this an even number
x_step=(x_max-x_min)/no_samples;

for n=1:no_samples+1
   
    diameters(n)=x_min+(n-1)*x_step;  
   %f(n)=1/sqrt(2*pi*sigma^2)*exp(-(x(n)-diameter)^2/(2*sigma^2));
    f(n)=exp(-(diameters(n)-diameter)^2/(2*sigma^2));
    
end

%normalized:
f2=f/sum(f);%so these are the weights for all the diameters.. 
%x's are the diameters

else
diameters=diameter;
f2=1;
end


%figure(1)
%plot(diameters,f2)

%%

%for, need to loop over sphere size:

%tic

i_parallel_cv=0;
i_perp_cv=0;
i_unpolarized_cv=0;


%for diameter_n=1:no_samples+1
for diameter_n=1:length(diameters) 
    

d=diameters(diameter_n);
a=d/2;

%from Bohren, rougly x number of terms (n) is sufficient.. in practice, let's use 100 terms, about 10x more 
%than we need, why not.

%silica in PDMS
%x=2*pi*N_pdms*a/lambda
%m=N_silica/N_pdms
%title_string_array={'Mie Scattering, 1.5\mum silica sphere in PDMS', '40x objectives, 575nm light'};

x=2*pi*n_0*a/lambda;
m=n_1/n_0;

%N_mix=1.4311;
%N_silica=1.4589;
%silica in glycerine
%{
x=2*pi*N_mix*a/lambda;
m=N_silica/N_mix;
%}

%title_text=strcat('Mie Scattering: ', num2str(diameter), ' um diameter ',sphere_n, ' in ', bulk_n, ', lambda=', num2str(lambda),' microns');   
title_text={['Scattering: ' num2str(diameter) ' um diameter ' sphere_n ' sphere in ' bulk_n],...
    ['$\lambda=$' num2str(lambda) ' um, ' pol_to_plot ' polarization, $n_0=$' num2str(n_0) ', $n_1=$' num2str(n_1)]};   

%how many terms, in Bessel expansion for Mie scattering, n
n_max=100;
S_1=0;
S_2=0;

%makes use of:
%function [ a_n,b_n ] = a_n_b_n_fn(m,x,n_max)
%function [pis ] = pi_tau_fn(theta,n_max )

%note we could make a look-up table for a_n and b_n 
a_n_vector=zeros(n_max,1);
b_n_vector=zeros(n_max,1);
Q_scat=0;

for n=1:n_max
    [a_n,b_n]=a_n_b_n_fn(m,x,n);   
    a_n_vector(n,1)=a_n;
    b_n_vector(n,1)=b_n; 
    %from Matzler document, this gives the scattering efficiency, angle
    %dependence already incorporated evidentally? I would like a little
    %more rigor here.
    Q_scat=Q_scat+2/x^2*(2*n+1)*(a_n*conj(a_n)+b_n*conj(b_n));
end

m=1:n_max;
scaling_vector=(2.*m+1)./(m.*(m+1));

%loop over theta:
angle_increment=.1;%deg
%angle_increment=1;
no_thetas=180/angle_increment;
delta_theta=angle_increment*pi/180;

%data
thetas=zeros(no_thetas,1);
S_vector=zeros(no_thetas,2);

counter=1;
for theta_deg=0:no_thetas
    theta=theta_deg/no_thetas*pi;
    thetas(counter,1)=theta;
    
    %calculate all the pi and taus one time through:
    pi_tau_list=pi_tau_fn(theta,n_max);
    pi_vector=pi_tau_list(:,1);
    tau_vector=pi_tau_list(:,2);
    
    %ok i think we should be able to do this as a dot product:
    S_vector(counter,1)=dot(scaling_vector'.*a_n_vector,pi_vector)+dot(scaling_vector'.*b_n_vector,tau_vector);
    S_vector(counter,2)=dot(scaling_vector'.*a_n_vector,tau_vector)+dot(scaling_vector'.*b_n_vector,pi_vector);
 
    %also let's generate the cosines and sines for later
    cosines(counter,1)=cos(theta);
    sines(counter,1)=sin(theta);

    counter=counter+1;
    
end

i_parallel=S_vector(:,2).*conj(S_vector(:,2));
i_perp=S_vector(:,1).*conj(S_vector(:,1));
%so about unpolarized light...
i_unpolarized=1/2*(i_parallel+i_perp);

%quantities of interest:
Q_s_unpolarized=1/2*(i_parallel-i_perp); %wait does q have angular dependence
S_11=i_unpolarized;
S_12=Q_s_unpolarized;


%with cv:
i_parallel_cv=i_parallel_cv+i_parallel*f2(diameter_n);
i_perp_cv=i_perp_cv+i_perp*f2(diameter_n);
i_unpolarized_cv=i_unpolarized_cv+i_unpolarized*f2(diameter_n);


end %diameters loop over.. 


%how long to process different sphere sizes.. takes a few seconds. 
%toc 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%                       what does it look like? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%variety of ways to display, select via hard-coding.


%1)Log_10 polar plot:
%{
figure(1)
polar(-thetas,log10(i_unpolarized/(min(i_unpolarized))),'m')
hold on
polar(thetas,log10(i_unpolarized/(min(i_unpolarized))),'m')
%polar(thetas,log10(i_perp),'b')
%polar(thetas,log10(i_parallel),'r')
%polar(-thetas,log10(i_perp),'b')
%polar(-thetas,log10(i_parallel),'r')
title('Mie Scattering: 1.54 um PS spheres in water, log_{10} scale')
%}

%2) linear polar plot:
%{
figure(2)
polar(-thetas,(i_unpolarized),'m')
hold on
polar(thetas,(i_unpolarized),'m')
polar(thetas,(i_perp),'b')
polar(thetas,(i_parallel),'r')
polar(-thetas,(i_perp),'b')
polar(-thetas,(i_parallel),'r')
%}

%3)flat plot (log)
%

%for paper, let's include normalization. 
norm_factor=1/(dot(i_unpolarized_cv,sin(thetas))*(thetas(2)-thetas(1))*2*pi);
norm_factor2=1/(dot(i_unpolarized_cv*norm_factor,sin(thetas))*(thetas(2)-thetas(1))*2*pi);

figure(18)
if strcmp(pol_to_plot,'circ')
semilogy(thetas*180/pi,(i_unpolarized_cv)*norm_factor,'-k','linewidth',.5)
varargout{2}=i_unpolarized_cv;
elseif strcmp(pol_to_plot,'co')
semilogy(thetas*180/pi,(i_parallel_cv)*norm_factor,'-k','linewidth',.5)
varargout{2}=i_parallel_cv;
elseif strcmp(pol_to_plot,'perp')
semilogy(thetas*180/pi,(i_perp_cv)*norm_factor,'-k','linewidth',.5)
varargout{2}=i_perp_cv;
elseif strcmp(pol_to_plot,'all')
semilogy(thetas*180/pi,(i_unpolarized_cv)*norm_factor,'-k','linewidth',.5)
hold on
semilogy(thetas*180/pi,(i_perp_cv)*norm_factor,'-k','linewidth',.5)
semilogy(thetas*180/pi,(i_parallel_cv)*norm_factor,'-k','linewidth',.5)
semilogy(thetas*180/pi,(i_parallel_cv+i_perp_cv)/2*norm_factor,'-k','linewidth',.5)

varargout{2}=(i_parallel_cv+i_perp_cv)/2*norm_factor;

end
title(title_text,'Interpreter','latex')
xlabel('\theta')
ylabel('Normalized Intensity')

hold on

varargout{1}=thetas*180/pi;


%semilogy(thetas*180/pi,(i_perp),'b')
%semilogy(thetas*180/pi,(i_parallel),'r')
%hold off
%title({'log_{10} Mie Scattering',title_string_description})
%legend('unpolarized','perp polarization','parallel polarization')
%}

%4) linear flat plot
%{
figure(5)
plot(thetas*180/pi,(i_unpolarized)/20,'m')
hold on
%plot(thetas*180/pi,(i_perp),'b')
%plot(thetas*180/pi,(i_parallel),'r')
%hold off
title('Mie Scattering')
legend('unpolarized','perp polarization','parallel polarization')
%}

%looking backward (log):
%{
figure(23)
plot((-thetas+pi)*180/pi,log10(i_unpolarized),'m')
hold on
plot((-thetas+pi)*180/pi,log10(i_perp),'b')
plot((-thetas+pi)*180/pi,log10(i_parallel),'r')
hold off
title('log_{10} Mie Scattering')
legend('unpolarized','perp polarization','parallel polarization')
xlim([0 110])
%}

%looking backward
%{
%figure(24)
plot((-thetas+pi)*180/pi,(i_unpolarized),'m')
hold on
plot((-thetas+pi)*180/pi,(i_perp),'b')
plot((-thetas+pi)*180/pi,(i_parallel),'r')
%hold off
title('Mie Scattering backward')
%legend('unpolarized','perp polarization','parallel polarization')
xlim([0 130])
%}

%{




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Finite point in back focal plane effects: light in collimated space is 
%at different angles, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%ok, so let's add a section that accounts for finite spot size in the back focus
%of the objective. We should be able to accomplish this by just creating
%shifted weighted versions and finding the sum..
%1)finite excitation area, not a problem, as the angles all map to the same 
%locations in the focal plane

%ok this is a little bit wrong, as we're just looking at a slice, we're not
%looking at the scattering pattern in the round. 

%ok first how big is the spot size at the input objective BFP? 
%85mm input lens, approx 5mm beam..
%w=lambda*f/(pi*w_0)

waist=525E-9*85E-3/(pi*5E-3)
%this gives a very small spot, ~3um

%need the objective focal lengths.. Mag=FL_tubelens/FL_obj
%Olympus FL_tubelens=180mm, 60x
%Leica FL_tubelens=200mm, 63x

%so 
FL_Olympus=180E-3/60
FL_Leica=200E-3/63

delta_theta_fwhm=5; %this angle spread can be related to spot size, delta_theta=delta_x/f
delta_theta_fwhm=0;
%delta_theta_fwhm=.6;%this corresponds to 100um spot size w/9mm fl (20x)input objective..(in deg)
%delta_theta_fwhm=5; %this is something like an 800um spot in the back focal plane
delta_theta_fwhm=waist/FL_Leica*180/pi%this is about .05 deg, 3um spot size on 60x objective..

delta_theta_fwhm=delta_theta_fwhm*100

sigma=delta_theta_fwhm/(2*sqrt(2*log(2)));

%sigma=0;

%construct a vector for 360 deg
thetas_360=[-flipud(thetas);thetas(2:end)];
scat_360=[flipud(S_11);S_11(2:end)];

if sigma~=0

for p=1:101
    phi=(p-51)/50*3*sigma;
    phis(p,1)=phi;
    phi_dist(p,1)=exp(-(phi)^2/(2*sigma^2));
        
end

phi_dist_normalized=phi_dist/(sum(phi_dist));
%figure
%plot(phis,phi_dist_normalized)

%figure
%plot(thetas_360,scat_360)

%interp to find shifted versions..
%figure(13)
%vq = interp1(x,v,xq)
total=0;
for q=1:101

phi=phis(q,1);    
%shifted(:,q)=interp1(thetas_360,scat_360,thetas_360+phi);
slice=interp1(thetas_360,scat_360,thetas_360+phi/180*pi)*phi_dist_normalized(q,1);
total=total+slice;
%plot(thetas_360,total)
%hold on


    
end

else
   total=scat_360;
end
%{
figure(6)
%thetas_36_prime=asin(N_medium*sin(observed_thetas));
semilogy(thetas_360*180/pi,total/max(total),'--b')
hold on
xlabel('angle in deg')
ylabel('normalized scattering intensity')
title('Intensity of scattering as fn of ouput angle for unpolarized input')

%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. Multiple scattering Estimate...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5-10-2016 work this over again..
%6-15-2016 and again... seems like normalization is off. 

%this will be in the round: 

%approach, let's do first very simple: use a combination of 1xscattered+2xscattered..

%one trick, how much will be scattered again? 
    %transmission_got=exp(-mu_scat_got*170)
    
scat_length=170;
transmission_got=exp(-1/scat_length*170);
scattered_fraction=1-transmission_got;

%simplest estimate, just assume the same fraction is scattered again. This
%overestimates b/c for those have already scattered, they have on average
%less than the full 170um thickness left, however, if they are scattered at
%angle, the effective thickness and opportunity to scatter again increased.
%Let's assume those issues are a wash and the same fracton scatters again.

    %scattered_fraction=.25;
scattered_again_fraction=scattered_fraction^2

%is it the relative fractions that matter? 

relative_scattered=scattered_fraction/(scattered_fraction+scattered_again_fraction)
relative_scattered_again=scattered_again_fraction/(scattered_fraction+scattered_again_fraction)
%because that's almost a quarter...


%%

thetas_to_blur=thetas*180/pi;
observed_intensity=i_unpolarized;%should I normalize this to make it a scattering fn?


th = (0:1:360);% this will be confusing as in our flattening scheme: theta is phi, 1 deg increments
r= (0:1:180); % r is theta. also 1 deg increments..increments are equally spaced..

%points we're going to look at in polar coordinates:
[TH,R] = meshgrid(th,r);
%translate to cartesian:
[X,Y] = pol2cart(TH,R);

%this finds the values by interpolating the x and y points in our mie
%radial plot, since it's cylindrically symmetrical, where r=sqrt(x^2+y^2)

%Z=interp1(theta angle, observed scattering intensity at that theta, length
%of R is theta for given angle- we don't care about the phi as it is azimuthally symetric)
Z=interp1(thetas_to_blur,observed_intensity/5.3104E7,sqrt(X.^2+Y.^2));

%for some reason there are some NaN's, kill them?:
Z(isnan(Z))=0;

%normalize so the area of the scattering distribution is 1:
%Z_normalized=Z/sum(sum(Z));
Z_normalized=Z;

%plot in 3D
%so this shows the forward scattering distribution stretched/mapped from a
%point to a plane
%{
figure(13)
%surf(X,Y,Z_normalized)
%surf(X,Y,log10(Z_normalized))
surf(X,Y,log10(Z))

shading flat
colormap summer
%}

%%
%wait a sec...-clean this up..
convolution_holder=zeros(size(Z));

%can we shift this to reflect this distribution but at a non-zero incident
%angle...?
%{
th_0=10*pi/180;
r_0=10;

x_0=r_0*cos(th_0);
y_0=r_0*sin(th_0);
Z2=interp1(forward_observed_thetas*180/pi,normalized_observed_intensity,sqrt((X-x_0).^2+(Y-y_0).^2));
%for some reason there are some NaN's, kill them:
Z2(isnan(Z2))=0;

%normalize so the area of the scattering distribution is 1:
Z2_normalized=Z2/sum(sum(Z));


%plot in 3D
%so this shows the forward scattering distribution stretched/mapped from a
%point to a plane
figure(14)
surf(X,Y,Z2_normalized)
shading flat
colormap summer

%the total is less for Z2 because some is lost beyond the NA that we can
%see. comparing:
%sum(sum(Z_normalized))
%sum(sum(Z2_normalized))

%yes, ok this puts the scattering distribution at another angle (rho and
%phi)

%}
%
%ok, so next, let's loop over our theta and phi angle grid...

%this is what we're doing, but we have them stored..
%{
for n_th=(1:1:360)*pi/180
   for n_phi=(1:1:90) 
   end
end
%}
total_for_normalization=sum(sum(Z)); %does this actually depend on the centering? maybe that's the problem..

%wrapped versions:
thetas_to_blur_extended=[thetas_to_blur;thetas_to_blur(2:end)+180];
observed_intensity_extended=[observed_intensity;flip(observed_intensity(1:end-1))]/5.3104E7;
        

%ok lets loop over all the angles in the grid..
for th_n=1:length(th)
    for r_n=1:length(r)
        
        %need the weight, i.e. fn(theta,phi) 
        weight=Z_normalized(r_n,th_n);
        
        %find the actual theta, phi, x and y for interp
        r_0=r(r_n);
        th_0=th(th_n);
        x_0=r_0*cos(th_0);
        y_0=r_0*sin(th_0);
        
        %this spits out the interpolated scattering pattern but shifted so
        %that its starting from the initial angle (at this point in the loop) in the angle grid..
            % Z2=interp1(thetas_to_blur,observed_intensity,sqrt((X-x_0).^2+(Y-y_0).^2));
     
            %%%%%%%%%%%%%%%%%%%%%%%%%hey do we need to
        %%%%%%%%%%%%%%%%%%%%%%%%%rap...%%%%%%%%%%%%%%
        % thetas to blur only go out to 180 deg, but could be larger no? 
        %let's make a thing that goes past 180 deg? actually the problem is
        %perhaps exagerrated looking at it on log scale:/.. wrapping is
        %more accurate though..
                 Z2=interp1(thetas_to_blur_extended,observed_intensity_extended,sqrt((X-x_0).^2+(Y-y_0).^2));

            %hmm does change the total but only very slightly, why is it
            %less than one now?
        
        %for some reason there are some NaN's, kill them:
        Z2(isnan(Z2))=0;
        %normalize so the area of the scattering distribution is 1:
        %so it won't be one, but less, because some is scattered beyond
        %NA.. for the shifted ones.
    %     Z2_normalized=Z2/total_for_normalization;
            Z2_normalized=Z2;
        %this adds the contribution from each tot he holder..weighted
        %by how much there would be scattered at that angle initially. 
        convolution_holder=convolution_holder+Z2_normalized*weight;
        
        
        
        
    end
    
    
end

%sum(sum(convolution_holder))


%now let's plot this business- the convolved version... 
%{
figure(15)
%surf(X,Y,convolution_holder)
surf(X,Y,log10(convolution_holder))
shading flat
colormap summer
%}
%ok great, it looks like a blurred version..


%%

%now we need relative fractions of each.
%seems like there may be something wrong with the scaling, as if the 2x
%scattered overall level is too large...




scaling_factor=1.713E4;
%scaling_factor=8347;

%figure(7)
%plot(r',Z_normalized(:,1)/max(Z_normalized(:,1)))
%semilogy(r',Z_normalized(:,1)/max(Z_normalized(:,1))*scaling_factor)

%hold on
%plot(r',convolution_holder(:,1)/max(convolution_holder(:,1)))
%semilogy(r',convolution_holder(:,1)/max(convolution_holder(:,1))*scaling_factor)

%xlim([0 50])


for relative_scattered_again=0:.2:1
Z_mix=relative_scattered_again*convolution_holder+Z_normalized*(1-relative_scattered_again);

%{
%plot of convololution in 3d
figure(16)
surf(X,Y,Z_mix)
shading flat
colormap summer
%}
figure(18)
%plot(r',Z_mix(:,1)/max(Z_mix(:,1)),'-k')
%semilogy(r',Z_mix(:,1)/max(Z_mix(:,1))*scaling_factor,'-k')
semilogy(r',Z_mix(:,1)/max(Z_normalized(:,1))*scaling_factor,'-k')
hold on
xlim([0 180])

end


%{
for relative_scattered_again=0:.25:1
Z_mix=relative_scattered_again*convolution_holder+Z_normalized*(1-relative_scattered_again);

%{
%plot of convololution in 3d
figure(16)
surf(X,Y,Z_mix)
shading flat
colormap summer
%}
figure(7)
%plot(r',Z_mix(:,1)/max(Z_mix(:,1)),'-k')
semilogy(r',Z_mix(:,1)/max(Z_mix(:,1))*scaling_factor,'-k')
xlim([0 180])

end
%}


%%

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of interest:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%see: http://omlc.org/education/ece532/class3/mie_exqs.html

%geometrical cross section
A=pi*a^2%microns^2

cross_sect=A*Q_scat%um^2
%scattering efficiency: 
%Qs/Ii=S_21: pg66?
concentration=.1;
percent_solution=.1;
volume_fraction=percent_solution/1.05%unitless
parts_h2o=50;
volume_fraction_adj=volume_fraction*1/(1+parts_h2o)

no_density=volume_fraction_adj/(4/3*pi*a^3);%/um^3
volume_frac=4/3*pi*a^3*concentration

%volume_frac=.001;
mu_scat=3/4*volume_frac/a*Q_scat
mu_scat=cross_sect*no_density%/um
%confusing about Q_scat, B&H discuss how to calculate,
%integrating/incorporating the angular dependence..
scat_length=1/mu_scat%um

    
%the anisotropy:
g=sum(S_11.*cosines.*sines*2*pi*delta_theta)/sum(S_11.*sines*2*pi*delta_theta)

%%

%calculate the dilution

%Q_scat

%required volume fraction
%fv=4/3*a/Q_scat*1/(1-g)*1/170% reduced scattering coefficient

fv=4/3*a/Q_scat*1/170%required volume fraction

%starting volume fraction given by ratio of masses: 
fv_initial=.1017/1.0525

v_in=1;

v_add=v_in*(fv_initial/fv-1)

%initial sample
fv_got=.0966/(1+3.95)

%next sample
fv_got=.0966/(1+80)
fv_got=.0966/(1+25)

%mu_scat=3/4*fv_got/a*Q_scat*(1-g)%reduced

mu_scat_want=3/4*fv/a*Q_scat%1/um
mu_scat_got=3/4*fv_got/a*Q_scat%1/um

%want:

scat_length_want=1/mu_scat_want%um
scat_length_got=1/mu_scat_got%um


transmission_got=exp(-mu_scat_got*170)


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of interest, take 2:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%see: http://omlc.org/education/ece532/class3/mie_exqs.html

%geometrical cross section
A=pi*a^2%microns^2, where a is the particle radius given above..

%note density of PS is 1.05g/cc, vs 1g/cc for H2O
cross_sect=A*Q_scat%um^2
%scattering efficiency: 
%Qs/Ii=S_21: pg66?

percent_solution=.1;%percent by mass
percent_solution=.097;
volume_fraction_0=percent_solution/1.05%unitless
parts_h2o=75;%good for .94um spheres
%parts_h2o=25;%good for 2um spheres
parts_h2o=50;

volume_fraction_adj=volume_fraction_0*1/(1+parts_h2o)%adjusted by adding H2O

no_density=volume_fraction_adj/(4/3*pi*a^3);%/um^3, no spheres per vol
    
mu_scat=cross_sect*no_density%/um

mu_scat_cm=mu_scat*10^4

%confusing about Q_scat, B&H discuss how to calculate,
%integrating/incorporating the angular dependence..
scat_length_um=1/mu_scat%um
    
%the anisotropy:
g=sum(S_11.*cosines.*sines*2*pi*delta_theta)/sum(S_11.*sines*2*pi*delta_theta)

%}

end





