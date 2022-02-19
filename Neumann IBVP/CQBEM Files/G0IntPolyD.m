function [ val ] = G0IntPolyD(a,b,s)
if s==b
valb=0; 
vala=(a-s).*(log(abs(a-s))-1);
elseif s==a
valb=(b-s).*(log(abs(b-s))-1);
vala=0;
else 
valb=(b-s)*(log(abs(b-s))-1);
vala=(a-s)*(log(abs(a-s))-1);  
end

val=valb-vala;
 
end

