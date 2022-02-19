function [MinEdge,MaxEdge,nEdgeHit]=FindRange(j,jv,LEdge,nVert,ThetaLoc,xva,yva,zva,xvb,yvb,zvb,Omega,xnodes,ynodes,znodes)

MaxEdge=ones(size(ThetaLoc))*(nVert+(jv-1));%;  
MinEdge=ones(size(ThetaLoc))*(jv+1);%
        
       
        for tv=1:nVert-2  %tv loops over all vertices can transmit to (excludes jv and jv+1), given by rv
             if jv<nVert
                 if tv<jv
                     rv=tv;   %when jv<nVert, transmit to 1:jv-1 and jv+2:nVert
                 else
                     rv=tv+2;
                 end
             else
                rv=tv+1;   %when jv=nVert, transmit to 2:nVert-1
             end
             
             %Gives directed angle between rays from ajj and bjj hitting each vertex and surface
             %normal using dot prod fmla and sin(pi/2 + q)=cos(q)
                    
            Anglea(tv)=asin(((xnodes(Omega(j,rv))-xva)*(xvb-xva)+(ynodes(Omega(j,rv))-yva)*(yvb-yva)+(znodes(Omega(j,rv))-zva)*(zvb-zva))/(sqrt((xnodes(Omega(j,rv))-xva)*(xnodes(Omega(j,rv))-xva)+(ynodes(Omega(j,rv))-yva)*(ynodes(Omega(j,rv))-yva)+(znodes(Omega(j,rv))-zva)*(znodes(Omega(j,rv))-zva))*LEdge));
            Angleb(tv)=asin(((xnodes(Omega(j,rv))-xvb)*(xvb-xva)+(ynodes(Omega(j,rv))-yvb)*(yvb-yva)+(znodes(Omega(j,rv))-zvb)*(zvb-zva))/(sqrt((xnodes(Omega(j,rv))-xvb)*(xnodes(Omega(j,rv))-xvb)+(ynodes(Omega(j,rv))-yvb)*(ynodes(Omega(j,rv))-yvb)+(znodes(Omega(j,rv))-zvb)*(znodes(Omega(j,rv))-zvb))*LEdge));
         
            MaxEdge=MaxEdge-(ThetaLoc>=Anglea(tv));%;  
            MinEdge=MinEdge+(ThetaLoc<=Angleb(tv));%
        end
        nEdgeHit=1+MaxEdge-MinEdge;
        
        MaxEdge=mod(MaxEdge,nVert);
        MinEdge=mod(MinEdge,nVert);
        
        MaxEdge(MaxEdge==0)=nVert;
        MinEdge(MinEdge==0)=nVert;
                
end 
       