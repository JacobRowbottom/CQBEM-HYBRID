function [xij yij Sbound] = findsplit(Lxj,slopeij,slopejj,xaj,yaj,xv,yv,ajj)
%
%xij are on the sending edge j
%xij, yij give the coordinates of the integration ranges within element je- in the simplest case they're just the element end points, if the end point of the ray changes receiving edge, they include the corresponding break points on the sending edge.
if abs(Lxj)==0
    %coords of intersection between ray and
    %integration element jj
    
    xij=xaj;                           
    yij=slopeij*(xij-xv)+yv;
                                                           
elseif abs(slopeij)>1e14
                              
    xij=xv;                        
    yij=slopejj*(xij-xaj)+yaj;
                                                           
else
    
    xij=(yaj-yv+slopeij*xv-slopejj*xaj)/(slopeij-slopejj);                            
    yij=slopejj*(xij-xaj)+yaj; 
                                                              
end

Sbound=ajj+sqrt((xij-xaj)^2+(yij-yaj)^2);

end

                           
