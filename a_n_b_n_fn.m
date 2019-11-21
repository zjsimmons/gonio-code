function [ a_n,b_n ] = a_n_b_n_fn(m,x,n_max)
%a_n_calc is a function to calculate the scattering coefficients a_n and
%b_n

%{
%order
n=1;

%relative index
N_1=1.1;
N=1;
m=N_1/N;

%wavelength
lambda=1E-6;
%sphere size
a=.5E-6;
%scaling parameter
x=2*pi*N*a/lambda;

%for water droplets, a la Bohren and Huffman:
x=3;
m=1.33+1i*10E-8;
%evaluating a_n's and b_n's for these parameters gives us what we'd expect.


%}

%a number of issues, first off, we need spherical bessel functions which
%are not built in, but we can relate to them through: 
% spherical bessel (n,z) = besselj(n+1/2,z)*sqrt(pi/(2*z)). We also need to
% generate derivative terms of the Ricatti functions. See notes. 

%let's try and get these going as anonymous functions of the form:
%sqr = @(x) x.^2;

psi =@(rho,n) rho*sqrt(pi/(2.*rho))*besselj(n+1/2,rho);

xi =@(rho,n) rho*sqrt(pi/(2.*rho))*(besselj(n+1/2,rho)+1i*bessely(n+1/2,rho));

psi_prime =@(rho,n) sqrt(pi/(2.*rho))*(rho.*besselj(n-1/2,rho)-n.*besselj(n+1/2,rho));

xi_prime =@(rho,n) sqrt(pi/(2.*rho))*(rho.*(besselj(n-1/2,rho)+1i*bessely(n-1/2,rho))...
    -n.*(besselj(n+1/2,rho)+1i*bessely(n+1/2,rho)));


%so we can put these all together to come up with the scattering
%coefficients, a_n and b_n 
%
n=n_max;

a_n=(m*psi(m*x,n)*psi_prime(x,n)-psi(x,n)*psi_prime(m*x,n))/...
    (m*psi(m*x,n)*xi_prime(x,n)-xi(x,n)*psi_prime(m*x,n));

b_n=(psi(m*x,n)*psi_prime(x,n)-m*psi(x,n)*psi_prime(m*x,n))/...
    (psi(m*x,n)*xi_prime(x,n)-m*xi(x,n)*psi_prime(m*x,n));
%}

%or make these anonymous too? 
%{
a_n_fn =@(m,x,n) (m*psi(m*x,n)*psi_prime(x,n)-psi(x,n)*psi_prime(m*x,n))/...
    (m*psi(m*x,n)*xi_prime(x,n)-xi(x,n)*psi_prime(m*x,n));

b_n_fn =@(m,x,n) (psi(m*x,n)*psi_prime(x,n)-m*psi(x,n)*psi_prime(m*x,n))/...
    (psi(m*x,n)*xi_prime(x,n)-m*xi(x,n)*psi_prime(m*x,n));
%}


%test
%{
for l=1:10
    as(l,1)=a_n_fn(m,x,l);
    bs(l,1)=b_n_fn(m,x,l);
    
end

as
bs
%}

end

