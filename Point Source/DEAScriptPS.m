
% CALCULATE THE COEFFICENTS GAMMA
stage =3 
Maxutilde1 = utildeB(:,pos); %Boundary solution at kfreq.
Maxutilde2 = CQBEM_par_Calc(k(pos+1),Ftilde_BC(:,pos+1),M,NVert,nEdge,cEdge,CPi,xi,yi,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny); %Boundary solution at wavenumber k(pos+1).

stage =3 
GammaScript; 

%CALCULATE THE BOUNDARY SOLUTION FOR HIGH FREQUENCIES VIA DEA

%%
stage=4
%Calculates Amplitudes via DEA
AmpS=zeros(length(VecInf2),N);
parfor ll=1:(N/2)+1  
     if abs(real(k(ll))) > abs(real(kfreq))%Freq 
        VecInf= DEA_Calc_PS(dx,k(ll),xv,yv,P0tilde(ll),PSx,PSy,NQ,Omega,MinDamp); %VecInf: boundary Density, Tj: No. of directions pointing into the domain along edge j                         
        Ampp = sqrt(abs(VecInf(ExtEdge>0)));        
        AmpS(:,ll) =Ampp;
    end
end
AmpS(:,N:-1:(N/2)+2)=AmpS(:,2:(N/2));

 s=zeros(1,M);     
 for ll=1:N 
      if abs(real(k(ll))) > abs(real(kfreq))%Freq     
         for i=1:NVert           
             s(CPGlobal(i)+1:CPGlobal(i+1))=CPi(CPGlobal(i)+1:CPGlobal(i+1)) - a(CPGlobal(i)+1);
             for iv=1:nEdge(i)
                 it=iv+cEdge(i);
                 utildeDEA(it,ll) =0;               
                 for n=1:NewTj(it) % Sums solutions of utildeDEA from each direction asocciated at coll pt it
                    utildeDEA(it,ll) = utildeDEA(it,ll) + AmpS(C(it)+n,ll)*exp(1i*(real(k(ll)))*(real(GammaS(i,n)) + sin(ThetaLocp(i,n))*(CPi(it)-CL(i))));
                 end               
             end
         end       
           utildeB(:,ll) = utildeDEA(:,ll);
     end
 end 