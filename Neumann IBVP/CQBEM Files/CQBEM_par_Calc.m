function [phihat, VRHS]=CQBEM_par_Calc(k,Ftilde_BC,M,NVert,nEdge,cEdge,CPi,xi,yi,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny)

   VSLP = zeros(M,M);  
   KDLP = zeros(M,M);

     for i=1:NVert
         for iv=1:nEdge(i)
             it=iv+cEdge(i);
             for j=1:NVert
                 for jv=1:nEdge(j)
                     jt=jv+cEdge(j);
                        if it==jt    
                            VSLP(it,jt)=(-1i/4)*(integral(@(q)VSLPIntPolyDiagonal(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt)) + (2i/pi)*G0IntPolyD(a(jt),b(jt),CPi(it)));
%                           KDLP(it,jt) = 0;
                        else
                            VSLP(it,jt)=(-1i/4)*integral(@(q)VSLPIntPoly(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j)),a(jt),b(jt));  
                            KDLP(it,jt)=(-1i/4)*integral(@(q)KIntPoly(q,k,xi(it),yi(it),xv(j),yv(j),CosEdgeAngle(j),SinEdgeAngle(j),CL(j),nx(jt),ny(jt)),a(jt),b(jt));      
                        end
                 end
             end
         end
     end
     KLHS = ((0.5)*eye(M) + KDLP); %LHS of transformed Helm eqs 
     VRHS = VSLP*Ftilde_BC; %RHS of transformed Helm eqs
     
 %%% - Rows contain coefficents of each wave number l.     
      phihat=KLHS\VRHS; 
 end