
% CALCULATE THE COEFFICENTS GAMMA
stage =3 
Maxutilde1 = utildeB(:,pos); %Boundary solution at kfreq.
Maxutilde2 = CQBEM_par_Calc(k(pos+1),Ftilde_BC(:,pos+1),M,NVert,nEdge,cEdge,CPi,xi,yi,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny); %Boundary solution at wavenumber k(pos+1).

GammaScriptL; 

%%CALCULATE THE BOUNDARY SOLUTION FOR HIGH FREQUENCIES VIA DEA
stage=4
%Calculates Amplitudes via DEA
AmpS=zeros(length(Amp1Final),N);
for ll=1:(N/2)+1  
    if abs(real(k(ll))) > abs(real(kfreq))
        VecInf =DEA_par_Calc_L(dx,k(ll),xv,yv,Ftilde_BC(:,ll),NQ,tht0,Omega,M,cEdge);
        Ampp = sqrt(VecInf(ExtEdge>0));        
%        Reorders Amplitudes for multi-domain DEA case
         Amp1 =Ampp(1:Tj(1,1)*nSideEls(1,1)); %edge 1    
         Amp3 =Ampp(Tj(1,1)*nSideEls(1,1)+1:sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2));
         Amp2 =Ampp(sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2)+1:end);
         Amp=[Amp1;Amp2;Amp3];

        AmpS(:,ll) = Amp;

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
                     utildeDEA(it,ll) = utildeDEA(it,ll) + AmpS(C(it)+n,ll)*exp(1i*(real(k(ll)))*((real(GammaS(i,n))) + sin(ThetaLocp(i,n))*CPi(it)));
                 end
             end
         end        
        utildeB(:,ll) = utildeDEA(:,ll); % Inputing the boundary solns from DEA into the overall boundary solution
    end
end 