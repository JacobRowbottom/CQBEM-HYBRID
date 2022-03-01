function [f, df]=WaveEqSol(x,alpha,t0)
     f=exp(-alpha*(x+t0).^2);
     df=-2*alpha*(x+t0).*exp(-alpha*(x+t0).^2);
end