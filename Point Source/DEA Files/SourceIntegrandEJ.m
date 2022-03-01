function [f] = SourceIntegrandEJ(s,xa,ya,za,xb,yb,zb,Ox,Oy,Oz,Damp,ThetaLoc,ThetaStep)

% Source Integrand- reflecting outer edge case 
% s from 0 to Dx

Dx=sqrt((xb-xa).^2+(yb-ya).^2+(zb-za).^2);

xs=xa+s.*(xb-xa)./Dx;
ys=ya+s.*(yb-ya)./Dx;
zs=za+s.*(zb-za)./Dx;

R0=sqrt((Ox-xs).^2+(Oy-ys).^2+(Oz-zs).^2);

SinTheta0=-((xb-xa).*(Ox-xs)+(yb-ya).*(Oy-ys)+(zb-za).*(Oz-zs))./(R0.*Dx);
CosTheta0=sqrt(1-SinTheta0.^2);
Theta0=asin(SinTheta0);
%pwc
CosThetaCut=(abs(ThetaLoc-Theta0)<0.5*ThetaStep)./CosTheta0;%
%pwl
%CosThetaCut=(abs(ThetaLoc-Theta0)<ThetaStep).*(ThetaStep-abs(ThetaLoc-Theta0))./(ThetaStep*CosTheta0);%

f=CosThetaCut./(sqrt(Dx));

end

