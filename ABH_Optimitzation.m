function [R,f]=ABH_Optimitzation(rvec,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the reflection coefficient of an Acoustic Black Hole (ABH)
% termination using transfer matrices that relate the acoustic pressure (P)
% and the acoustic velocity (Q)
%
% Inputs:
%   rvec = vector to stablish the inner radius of the black hole (in tan per cert with respect to R) for each ring.
%
% Outputs:
%   R = reflection coefficient
%   f = frequency vector
%
% Fixed variables:
%   fm = sampling frequency
%   N = number of rings
%   hring = ring thickness
%   xl = distance from the end termination for the truncation (m)
%   damping = damping value, typically 0.05
%   geo = 'a1' --> used to change the geometry. Typically set to the case named 'a1'
%   conf = 'lin','qua','vec' --> linear ABH, quadratic ABH or with the profile introduced by a vector
%
% Examples:
% >> ABH_Optimitzation(20:1.5:80);
%
% Author: Marc Arnela, 2021
% Last update: 20 May 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=40;
hring=0.001/1000;
xl=1.0e-3;


if strcmp(type,'vec')
    Var=ABH_ABCD_PQesq (1000,N,hring,xl,0.05,'a1','vec',rvec); 

else
    Var=ABH_ABCD_PQesq (1000,N,hring,xl,0.05,'a1','qua',rvec); 
end
%Var=ABH_ABCD_PQesq (1000,N,hring,xl,0.05,'a1','vec',rvec,'kz',100); 
%Var=ABH_ABCD_PQesq (1000,N,hring,xl,0.05,'a1','vec',rvec); 
%Var=ABH_ABCD_PQesq (1000,N,hring,xl,0.05,'a1','lin',rvec); 

R=Var.R;
f=Var.f;

if nargout==0
    %set(0,'DefaultTextInterpreter', 'tex');

    %fontsize = 14;
    %set(0,'defaultaxesfontsize',fontsize);
    %set(0,'defaulttextfontsize',fontsize);
    %figure;
    %ZinN=Var.ZinN;
    %subplot(211), plot(f/1000,abs(ZinN)); ylabel('$$|Z_{\rm in}|$$','Interpreter','latex');xlabel('Frequency (kHz)','Interpreter','latex'); ylim([0 10])
    
    %subplot(212), plot(f,abs(R)); ylabel('$$|\mathcal{R}|$$','Interpreter','latex');  xlabel('Frequency (kHz)','Interpreter','latex'); ylim([0 1])
    
end    