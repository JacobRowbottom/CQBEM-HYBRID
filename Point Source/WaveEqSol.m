function [f, df]=WaveEqSol(x,n,t0)
 
     f=1+cos(n*pi*x+pi);   
     df=-n*pi*sin(n*pi*x+pi);
     
%      f=exp(-n*(x+t0).^2);
%      df=-2*n*(x+t0).*exp(-n*(x+t0).^2);

end