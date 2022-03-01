function [DX]=FindGrid(j,nVert,LEdge,Lx,Ly,Lz,xRay,yRay,zRay,xva,yva,zva,xvb,yvb,zvb,xnodes,ynodes,znodes,Omega,ThetaLoc,MaxEdge,nEdgeHit)

        MinElts=max(nEdgeHit);

         %Need to divide up edge jv using
        %points where hit vertices
        nL=length(ThetaLoc);
        
        xSplit=NaN*ones(max(nEdgeHit)-1,nL); 
        ySplit=NaN*ones(max(nEdgeHit)-1,nL);        
        zSplit=NaN*ones(max(nEdgeHit)-1,nL);        
            
               
        for iSint=1:MinElts-1 
                        
         %   iVertpre=MaxEdge-iSint+1;
            
            iVertp=mod(MaxEdge-iSint+1,nVert);
            iVertp(iVertp==0)=nVert;
         %   iVertp(nEdgeHit==1)=nVert(j)+1;
            
            xv=xnodes(Omega(j,iVertp));
            yv=ynodes(Omega(j,iVertp));
            zv=znodes(Omega(j,iVertp));
            
            [xSplit(iSint,iSint<=nEdgeHit-1),ySplit(iSint,iSint<=nEdgeHit-1),zSplit(iSint,iSint<=nEdgeHit-1)]=Intersection3D(xai,yai,zai,xbi,ybi,zbi,xaj,yaj,zaj,xbj,ybj,zbj)
            
%             if Lx==0
%             %coords of intersection between ray and
%             %start edge
%                 xSplit(iSint,iSint<=nEdgeHit-1)=ones(size(xv(iSint<=nEdgeHit-1)))*xva;                           
%                 ySplit(iSint,iSint<=nEdgeHit-1)=slopeij(iSint<=nEdgeHit-1).*(xSplit(iSint,(iSint<=nEdgeHit-1))-xv(iSint<=nEdgeHit-1))+yv(iSint<=nEdgeHit-1);
%                
%                                                            
%             elseif abs(slopeij)>1e14
%                 xSplit(iSint,iSint<=nEdgeHit-1)=xv(iSint<=nEdgeHit-1);                        
%                 ySplit(iSint,iSint<=nEdgeHit-1)=slope*(xSplit(iSint,(iSint<=nEdgeHit-1))-xva)+yva;
%                 
%             else
%                 xSplit(iSint,iSint<=nEdgeHit-1)=(yva-yv(iSint<=nEdgeHit-1)+slopeij(iSint<=nEdgeHit-1).*xv(iSint<=nEdgeHit-1)-slope*xva)./(slopeij(iSint<=nEdgeHit-1)-slope)        ;                    
%                 ySplit(iSint,iSint<=nEdgeHit-1)=slope*(xSplit(iSint,(iSint<=nEdgeHit-1))-xva)+yva; 
%                 
%             end
    
        end     
       %      
      %% Should generalise the below as only works for up to 3 target edges, not an immediately major issue though...

        DX=zeros(MinElts,length(nEdgeHit));
        
        DX(1:MinElts,nEdgeHit==1)=LEdge/MinElts;
        %No splitting needed so just equispace the MinElts elements
        
        DX(1,nEdgeHit==2)=sqrt((xSplit(1,nEdgeHit==2)-xva).^2+(ySplit(1,nEdgeHit==2)-yva).^2);
        DX(MinElts,nEdgeHit==2)=sqrt((xvb-xSplit(1,nEdgeHit==2)).^2+(yvb-ySplit(1,nEdgeHit==2)).^2);
       
        %In case of one or more split take the 1st and last elements either
        %side of the split
        if MinElts>2
                    
            [MaxL,MaxP]=max(DX); % Refine the longest edge
                       
            DX(1,MaxP==1 & nEdgeHit==2)=0.5*MaxL(MaxP==1 & nEdgeHit==2);
            DX(2,nEdgeHit==2)=0.5*MaxL(nEdgeHit==2);
            DX(3,MaxP==3 & nEdgeHit==2)=0.5*MaxL(MaxP==3 & nEdgeHit==2);
       
        
            DX(1,nEdgeHit==3)=sqrt((xSplit(1,nEdgeHit==3)-xva).^2+(ySplit(1,nEdgeHit==3)-yva).^2);
            DX(2,nEdgeHit==3)=sqrt((xSplit(1,nEdgeHit==3)-xSplit(2,nEdgeHit==3)).^2+(ySplit(1,nEdgeHit==3)-ySplit(2,nEdgeHit==3)).^2);
            DX(3,nEdgeHit==3)=sqrt((xSplit(2,nEdgeHit==3)-xvb).^2+(ySplit(2,nEdgeHit==3)-yvb).^2);
        end


end

