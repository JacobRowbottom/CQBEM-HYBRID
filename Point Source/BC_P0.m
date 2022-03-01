P0t=zeros(N);
for nn=1:N
    l=(lambda^(nn-1));
    P0t(nn)=l*P0(t(nn),alpha,t0);
end
P0tilde=fft(P0t);

ri=sqrt((xi-PSx).^2+(yi-PSy).^2);

Bessel = (-k).*besselh(1,1,k.*ri);

Ker = (1./ri).*(nx.*Bessel.*(xi-PSx) +  ny.*Bessel.*(yi-PSy));

 for nn=1:N
     for i=1:NVert
         for iv=1:nEdge(i)
             it=iv+cEdge(i);
             Ftilde_BC(it,nn)=(-1i/4)*P0tilde(nn)*Ker(nn,it);
         end
     end
 end

