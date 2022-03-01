function [xTarg yTarg] = FindTarget(xStart,yStart,ThetaLoc,slopeij)

if abs(Lx(j,MinEdge))==0 %Spec. case: receiving edge is vertical

    xji(1+nEdgeHit)=xv(MinEdge);
    yji(1+nEdgeHit)=slopeij*(xji(1+nEdgeHit)-xbj)+ybj;
%                             
elseif abs(slopeij)>1e14 %Spec. case: ray is near vertical
    slopemin=Ly(MinEdge)/Lx(MinEdge);
    xji(1+nEdgeHit)=xbj;
    yji(1+nEdgeHit)=slopemin*(xji(1+nEdgeHit)-xv(MinEdge))+yv(MinEdge);         
%                             
else % general case                              
                  
    slopemin=Ly(MinEdge)/Lx(MinEdge);
    xji(1+nEdgeHit)=(yv(MinEdge)-ybj+slopeij*xbj-slopemin*xv(MinEdge))/(slopeij-slopemin);
    yji(1+nEdgeHit)=slopemin*(xji(1+nEdgeHit)-xv(MinEdge))+yv(MinEdge); 
%                                 
end
%                           %-------------------------------------------------
%Compute coords of intersection between ray from ajj on elt je and final edge MaxEdge 
%                           
if abs(Lx(MaxEdge))==0
    xji(1)=xv(MaxEdge);
    yji(1)=slopeij*(xji(1)-xaj)+yaj;
                          
elseif abs(slopeij)>1e14
    slopemax=Ly(MaxEdge)/Lx(MaxEdge);
    xji(1)=xaj;
    yji(1)=slopemax*(xji(1)-xv(MaxEdge))+yv(MaxEdge);                                   
        
else
    slopemax=Ly(MaxEdge)/Lx(MaxEdge);
    xji(1)=(yv(MaxEdge)-yaj+slopeij*xaj-slopemax*xv(MaxEdge))/(slopeij-slopemax);
    yji(1)=slopemax*(xji(1)-xv(MaxEdge))+yv(MaxEdge); 
end

%Need to divide up position integration range [ajj, bjj] at
%points where hit vertices - list of intergation
%boundary points is Sbound
                         
Sbound(1)=ajj;
Sbound(1+nEdgeHit)=bjj;
xij(1)=xaj;
xij(1+nEdgeHit)=xbj;
yij(1)=yaj;
yij(1+nEdgeHit)=ybj;
%Compute coords of intersection between ray and any intermediate vertices mapped to from elt je
%together with corresponding integration splitting points
          
for iSint=1:nEdgeHit-1 
    iVertp=MaxEdge-iSint+1;
    if iVertp<1
       iVertp=iVertp+nVert;
    end                               
    xji(iSint+1)=xv(iVertp);
    yji(iSint+1)=yv(iVertp);
end


end

