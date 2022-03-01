function [f] = P0(t,alpha,t0)

f=sqrt(alpha/pi).*sin(4*pi*t).*exp(-alpha*(t-t0).^2); %Given by equation (3.42) 
% f=sqrt(alpha/pi).*exp(-alpha*(t-t0).^2); % Equation (5.34) in thesis.

end

