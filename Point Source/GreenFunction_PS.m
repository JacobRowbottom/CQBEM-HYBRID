function [G]=GreenFunction_PS(IPx,IPy,PSx,PSy,c,tn,alpha,t0)       
             
rIPPS = sqrt((IPx - PSx)^2 + (IPy - PSy)^2); % distance from interior pts and PS. 
         
%Green function for homogeneous eqn for PS prob in the time domain.
G = 0.5*(sqrt(alpha)/(pi^(3/2)))*integral(@(tau)GreenFunWaveEQ(tn,tau,rIPPS,alpha,t0,c),0,max(0,tn-rIPPS/c)); 
