function [Int ] = KIntPoly(q,k,xpi,ypi,qxv,qyv,CosEdgeAngle,SinEdgeAngle,CL,nxq,nyq)

qx=qxv+(q-CL)*CosEdgeAngle; %cartesian coordinates in terms of arclength
qy=qyv+(q-CL)*SinEdgeAngle;

D=sqrt((xpi-qx).^2+(ypi-qy).^2);
Bessel = (-k).*besselh(1,1,k.*D);

Ker = (1./D).*(nxq.*Bessel.*(qx-xpi) +  nyq.*Bessel.*(qy-ypi));

Int=Ker;

end