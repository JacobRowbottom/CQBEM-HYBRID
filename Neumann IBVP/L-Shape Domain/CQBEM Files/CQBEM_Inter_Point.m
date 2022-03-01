function [vtildeInt]=CQBEM_Inter_Point(k,vtildeB,vBC,M,IPx,IPy,NVert,nEdge,cEdge,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nxq,nyq)

VSLPint = zeros(M,1);  
KDLPint = zeros(M,1);


for i=1:NVert
    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        VSLPint(it)=(-1i/4)*integral(@(q)VSLPIntPoly(q,k,IPx,IPy,xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i)),a(it),b(it));     
        KDLPint(it)=(-1i/4)*integral(@(q)KIntPoly(q,k,IPx,IPy,xv(i),yv(i),CosEdgeAngle(i),SinEdgeAngle(i),CL(i),nxq(it),nyq(it)),a(it),b(it));            
    end
end

vtildeInt = (vBC.')*VSLPint - (vtildeB.')*KDLPint; %This is correct formulation for when Gk = -i/4 coefficent.
