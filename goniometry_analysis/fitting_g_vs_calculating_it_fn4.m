function [ g_via_integration,g_via_fit,scattered_bw,r_squared,a,b ] = fitting_g_vs_calculating_it_fn4( angle_vector,theta_lower_limit,theta_upper_limit,scattering_vector,weights_in )

%DESCRIPTION:
%This function generates anisotropy values by both integration and by
%fitting HG. This most-recent, improved version will spit out fitting
%characteristics as well as takes angle inputs. These angle inputs set the
%range over which the HG fit is done, this is useful so we can easily see
%how the fit changes with excluding more or less very forward contribution,
%as well as so the same program can be used for partial as well as to 180
%degree traces.

%INPUTS:
%angle_vector- list of input angles
%theta_lower_limit - lower limit of angles to include in HG fit
%theta_upper_limit - upper limit of angles to include in HG fit
%scattering_vector - signal trace for those angles
%weights - is to do weighted least squares rather than just least squares

%note re: weights: better for our phase function data since the scattering
%is so forward directed. The weights correspond to avg
%1/var(different scans). Averaged together data from different colors to
%get a smoother curve. These weights are determined elsewhere from our
%experimental data. Entering 0 rather than a weights vector will result in
%no weighting being used.

%OUTPUTS:
%g_via_integration
%g_via_fit
%r_squared - goodness of fit
%a - fit number
%b - fit number

%zjs
%cleaned up 2016.12.9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if weights_in==0
    weights_in=ones(size(angle_vector));
end

bin_vals=scattering_vector;
bins=angle_vector;
angle_step_size=.1;

bin_vals(isnan(bin_vals))=0;

%re: calc g via integration
integrands = bsxfun(@times,bin_vals',sind(bins)')*2*pi*angle_step_size;
integrals=sum(integrands);
normalization_factors=1./integrals; %1xno bins row % norm factors for each phi bin
bin_vals_normalized=bsxfun(@times,bin_vals,normalization_factors');%this multiplies
%each all the rows - p(thetas) by the appropriate norm factor for that
%row (phi bin)

%figure(13)
%plot(bin_vals_normalized)
%size(bin_vals_normalized)
scattered_bw=sum(bin_vals_normalized(bins>=175));


%next calculate g's for each phi bin, multiply the p(thetas) by the
%appropriate stuff to get the integral for g's:
integrands2 = bsxfun(@times,bin_vals_normalized',(sind(bins).*cosd(bins))')*2*pi*angle_step_size;
g_via_integration=sum(integrands2);

%re: fit g
%{
theta_lower_limit=5;%brain
theta_lower_limit=10;
theta_lower_limit=4.5;
theta_upper_limit=181;
%theta_upper_limit=50;
%}

y=bin_vals_normalized(bins>theta_lower_limit&bins<theta_upper_limit);
weights=weights_in(bins>theta_lower_limit&bins<theta_upper_limit);
x=bins(bins>theta_lower_limit&bins<theta_upper_limit);

%figure(91)
%plot(bins,bin_vals_normalized)
%drawnow
%hold on

x_spliced=[-flip(x);x];
y_spliced=[flip(y);y];
weights_spliced=[flip(weights);weights];

%size(x_spliced);
%size(weights_spliced);

x=x_spliced;
y=y_spliced;
%y=log(y_spliced);

hg=fittype('1/(4*pi)*(1-g^2)/(1+g^2-2*g*cos(x*pi/180))^(3/2)*a+b',...
    'dependent',{'y'},'independent',{'x'},'coefficients',{'g','a','b'});

%ok, that didn't work well: 
%hg=fittype('1/(4*pi)*((1-g^2)/(1+g^2-2*g*cos(x*pi/180))^(3/2)*(1-b)+b)',...
%    'dependent',{'y'},'independent',{'x'},'coefficients',{'g','b'});


%hg=fittype('log(1/(4*pi)*(1-g^2)/(1+g^2-2*g*cos(x*pi/180))^(3/2)*a+b)',...
%     'dependent',{'y'},'independent',{'x'},'coefficients',{'g','a','b'})

%re WM function
%th is theta
%k is wavenumber
%lc is a lengthscale representing the 'outer' scale or where the fractal turns over to exponential.
%m is an old notation parameter.. mass fractal dimension D = 2m.

%th=x_spliced;

%wm=fittype('(k^6*lc^6*(1+4*k^2*lc^2)^m*(m-3)*(m-2)*(m-1)*(1-2*k^2*lc^2*(-1+cos(th))).^-m.*(3+cos(2*th)))/((1+4*k^2*lc^2)^m*(1+2*k^2*lc^2*(-1+2*k^2*lc^2*(m-2))*(m-3))*pi-(1+4*k^2*lc^2)*(1+2*k^2*lc^2*(1+m)+4*k^4*lc^4*(4+(m-3)*m))*pi)',...
% 'dependent',{'y'},'independent',{'th'},'coefficients',{'k','lc','m'})

%manually specified bounds:
options = fitoptions(hg);
%options.Lower = [.5 0.0125 0];
%options.Upper = [.99 0.0125 1];
%brain: grey matter
%options.Lower = [.5 0 0];
%options.Upper = [.99 1 1];
%retina:
%options.Lower = [.5 -1 -1];
%options.Upper = [.9999 1 1];
%white matter
options.Lower = [.5 0 0];
options.Upper = [.9999 10 0];

%options.Lower = [.5 -5 0];
%options.Upper = [.9999 5 0];

options.Weights=weights_spliced;

[myfit,gof]=fit(x,y,hg,options)

%[myfit,gof]=fit(th,y,wm,options)
r_squared=gof.rsquare;
    a=myfit.a;
b=myfit.b;
g_via_fit=myfit.g;

figure(17)
plot(myfit,x,y)
hold on

end
