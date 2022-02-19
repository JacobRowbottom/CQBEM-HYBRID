%% Calculates Interior solution via CQBEM using the Boundary solution previously calculated via the Hybrid method.
parfor ll =1:N 
    tn = (ll-1)*dt;
    for ipx=1:NP
         for ipy=1:NP  
             utildeInter(ipx,ipy,ll) = CQBEM_Inter_Point(k(ll),utildeB(:,ll),Ftilde_BC(:,ll),M,IPx(ipx),IPy(ipy),NVert,nEdge,cEdge,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny);
         end
    end
end

% INVERSE FFT FOR INTERIOR SOLUTION  
InputInt=zeros(NP,NP,N);
uInter=zeros(NP,NP,N);
for ipx=1:NP
    for ipy=1:NP
        InputInt(ipx,ipy,:) = (ifft(utildeInter(ipx,ipy,:))); 
    end
end
for n=1:N
    uInter(:,:,n)=((lambda)^(1-n))*InputInt(:,:,n); % uInter = Interior solution in time domain at each interior point
end

%% Exact Solution 
% for ipx=1:NP 
%     for ipy=1:NP
%         for kk=1:N
%             tn = (kk-1)*dt;  
%             uExactInt(ipx,ipy,kk)=0.5*(WaveEqSol(IPx(ipx)-c*tn,BCc,BCt0)+ WaveEqSol(IPx(ipx)+c*tn,BCc,BCt0));
%         end
%     end
% end

 for tt=1:N
     uInterior(tt) = real(uInter(1,1,tt));
%      uExact(tt) = real(uExactInt(1,1,tt));
 end

 figure
 plot(t,uInterior)
 
