function [f] = GreenFunWaveEQ(t,tau,r,alpha,t0,c)
 
f=zeros(size(tau)); 
f((t-tau-r/c)>0)=sin(4*pi*tau((t-tau-r/c)>0)).*exp(-alpha*(tau((t-tau-r/c)>0)-t0).^2)./sqrt((t-tau((t-tau-r/c)>0)).^2-(r/c)^2);
% f((t-tau-r/c)>0)=exp(-alpha*(tau((t-tau-r/c)>0)-t0).^2)./sqrt((t-tau((t-tau-r/c)>0)).^2-(r/c)^2);

end