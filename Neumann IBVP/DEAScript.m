
% CALCULATE THE COEFFICENTS GAMMA
stage =3 
Maxutilde1 = utildeB(:,pos); %Boundary solution at kfreq.
Maxutilde2 = CQBEM_par_Calc(k(pos+1),Ftilde_BC(:,pos+1),M,NVert,nEdge,cEdge,CPi,xi,yi,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny); %Boundary solution at wavenumber k(pos+1).

GammaScript; 

%CALCULATE THE BOUNDARY SOLUTION FOR HIGH FREQUENCIES VIA DEA

%%
stage=4
%Calculates Amplitudes via DEA
AmpS=zeros(length(VecInf2),N);
parfor ll=1:(N/2)+1  
    if abs(real(k(ll))) > abs(real(kfreq))
        VecInf =DEA_par_Calc(dx,k(ll),xv,yv,Ftilde_BC(:,ll),NQ,tht0,Omega);        
        Ampp = sqrt(VecInf(ExtEdge>0));       
        AmpS(:,ll) = Ampp;
    end
end
AmpS(:,N:-1:(N/2)+2)=AmpS(:,2:(N/2));

 
 %%Calculates Equation (5.21) in Thesis. 
for ll=1:N 
    if abs(real(k(ll))) > abs(real(kfreq))
         for i=1:NVert
             for iv=1:nEdge(i)
                 it=iv+cEdge(i);
                 utildeDEA(it,ll) =0;
                 for n=1:NewTj(it) %Solution for wave on left edge % Sums solutions of utildeDEA from each direction asocciated at coll pt it.
                     utildeDEA(it,ll) = utildeDEA(it,ll) + AmpS(C(it)+n,ll)*exp(1i*(real(k(ll)))*((real(GammaS(i,n))) + sin(ThetaLocp(i,n))*CPi(it)))*exp(-imag(GammaS(i,n))*real(k(pos+1))*(real(k(pos+1)))/real(k(ll)));
                 end
             end
         end      
        utildeB(:,ll) = utildeDEA(:,ll); % Inputing the boundary solns from DEA into the overall boundary solution
    end
end 