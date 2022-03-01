function [df]=F_BC(xi,yi,tht0,tn,c,alpha,t0)
% xi =0, since the boundary condition is only along the edge where x=0. 

[f, dfm]=WaveEqSol(xi*cos(tht0)+yi*sin(tht0)-c*tn,alpha,t0);
[f, dfp]=WaveEqSol(xi*cos(tht0)+yi*sin(tht0)+c*tn,alpha,t0);

df = (0.5)*dfm; %(0.5)*(dfm + dfp); 

end


