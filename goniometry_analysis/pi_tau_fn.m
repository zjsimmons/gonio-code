%function [pi_output,tau_output,pis ] = pi_fn(theta,n_max )
function [pis ] = pi_tau_fn(theta,n_max )

%pi_fn
%zjs 1.10.2015
%This function calculates the pi and tau legendre-function based functions
%necessary for the scattering problem. This version returns all the values
%up to n_max for a given theta in an array.

%small values
pi_0=0;
pi_1=1;
tau_1=cos(theta)*pi_1;

%vector of all values: column 1: pis, column 2: taus
pis=zeros(n_max,2);
pis(1,1)=pi_1;
pis(1,2)=tau_1;

%if you're looking for the n=1 term
if n_max==1
    pi_output=pi_1;
    tau_output=tau_1;
else
    
    %otherwise all the other terms are derived recursively
    pi_n_minus_1=pi_1;
    pi_n_minus_2=pi_0;
    
    for q=2:n_max
        pi_output=(2*q-1)/(q-1)*cos(theta)*pi_n_minus_1-q/(q-1)*pi_n_minus_2;
        tau_output=q*cos(theta)*pi_output-(q+1)*pi_n_minus_1;
        pi_n_minus_2=pi_n_minus_1;
        pi_n_minus_1=pi_output;        
        pis(q,1)=pi_output;
        pis(q,2)=tau_output;
    end
    
end



end

