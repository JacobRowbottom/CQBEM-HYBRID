function [Int ] = VSLPIntPoly(q,k,xpi,ypi,qxv,qyv,CosEdgeAngle,SinEdgeAngle,CL)

qx=qxv+(q-CL)*CosEdgeAngle; %cartesian coordinates in terms of arclength
qy=qyv+(q-CL)*SinEdgeAngle;

D=sqrt((xpi-qx).^2+(ypi-qy).^2);

Ker = besselh(0,1,k*D);

Int=Ker;

end