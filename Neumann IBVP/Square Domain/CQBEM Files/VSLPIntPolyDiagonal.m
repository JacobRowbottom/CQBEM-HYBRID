function [Int ] = VSLPIntPolyDiagonal(q,k,xpi,ypi,qxv,qyv,CosEdgeAngle,SinEdgeAngle,CL)
E=double(eulergamma);
IntGrandpre=ones(size(q)).*(1+(2i/pi)*(E+log(k/2)));

qx=qxv+(q-CL)*CosEdgeAngle; %cartesian coordinates in terms of arclength
qy=qyv+(q-CL)*SinEdgeAngle;

D=sqrt((xpi-qx).^2+(ypi-qy).^2);

Ker = besselh(0,1,k*D) - (2i/pi)*log(D);
Ker(D<1e-15)=IntGrandpre(D<1e-15); 

Int=Ker;

end