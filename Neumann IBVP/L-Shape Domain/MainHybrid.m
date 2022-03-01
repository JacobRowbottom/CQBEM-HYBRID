clear
tic

%% Prelimnaries and set up 
SL=zeros(1,6);
nEdge=zeros(1,6);
h=zeros(1,6);
CosEdgeAngle=zeros(1,6);
SinEdgeAngle=zeros(1,6);
CL=zeros(1,6);
cEdge=zeros(1,6);

%Parameters
DEA=1; %DEA =1 indicated DEA based hybrid method. DEA=0 indicaates SHFA Hybrid method.
N=64; %No. of time steps
M=16; %No. of Boundary elements 
NQ=2; %Define global direction set and elements multiples of 4 give symmetry in each quadrant
alpha=6^2; %Boundary condition constant that controls the bandwidth
t0=1; %Boundary condition constant that controls the peak of the wave
T = 2; %Time interval
Freq = 10; % Frequnecy switch value 
tht0=0; %Plane wave BC incoming angle - only works for angles from 0 up to strictly less than pi/2

c=1; %Wave speed
dt=T/N; %%% time step size 
t = linspace(t0,T,N); %Time Steps
lambda= 10^(-8/N); %dt^(3/(N-1)) % lambda: radius of contour 

NP=1; %number of interior points
IPx=0.5; %x -Coordinate
IPy=0.5; %y- Coordinate

%% Define vertices for domains and edge split for DEA

% L-shaped Domain
 xv = [0 1 1 0.5 0.5 0]; 
 yv = [0 0 1 1 0.5 0.5];
Omega=[1 2 5 6; 2 3 4 5]; %Defines the edge no. of each sub-domain for DEA

NVert=length(xv); %No. of vertices of domain 
for j=1:NVert
    jp = j+1;
    if (jp > NVert)
        jp=jp-NVert;
    end
    SL(j)=sqrt((xv(j)-xv(jp))^2+(yv(j)-yv(jp))^2); %Side Length of edges of domain
    CL(j)=sum(SL(1:j-1)); %Cummalative length of edges of domain
    CosEdgeAngle(j)=(xv(jp)-xv(j))/SL(j);
    SinEdgeAngle(j)=(yv(jp)-yv(j))/SL(j);
    
end

L=sum(SL); %Length/perimeter of polygon
dx = L/M; %Average Element size
x=dx:dx:L; 

for j=1:NVert
    nEdge(j) = round(SL(j)/dx); %Number of elements on each side
    cEdge(j)=sum(nEdge(1:j-1)); % cummaltive number of elements on edge
    h(j)=SL(j)/nEdge(j); %Stepsize on each side
    
   for je=1:nEdge(j) 
       % Calculating element endpoint coordinates. 
        a(je+cEdge(j))=CL(j)+(je-1)*h(j); 
        b(je+cEdge(j))=CL(j)+je*h(j);
       
        ax(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je-1);
        ay(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je-1); 

        bx(je+cEdge(j))=xv(j)+h(j)*CosEdgeAngle(j)*(je);
        by(je+cEdge(j))=yv(j)+h(j)*SinEdgeAngle(j)*(je);
              
   end
   %Calculating unit normal vectors
   for je=1:nEdge(j) 
       jt = je +cEdge(j);       
       r = sqrt((bx(jt)-ax(jt))^2 + (by(jt)-ay(jt))^2); 
       txq(jt) = (bx(jt) - ax(jt))/r;
       tyq(jt) = (by(jt) - ay(jt))/r;   
   end
end
% cartesian coordinates of unit normal vector of q.
nx = -tyq; 
ny = txq;  

%collocation points 
CPi= 0.5*(a+b); 
%collocation cartesian coordinates
xi= 0.5*(ax+bx); 
yi= 0.5*(ay+by); 

%% Calculate the complex wave numbers using an A-stable Mulitstep method.
zeta_l=zeros(N,1);
for m=1:N
    z=lambda*exp((-2i*pi*(m-1))/(N)); %calculates the index numbers  
    BDF2= (0.5)*(z^2 -4*z+3); 
%     BE = 1-z; %Backward Euler 
%     Trap = 2*(1-z)/(1+z); %Trapezium
    zeta_l(m)=BDF2/dt; %complex wave numbers from l=0:N-1. 
end
k=1i*zeta_l; %The complex wave numbers

addpath(genpath('CQBEM Files'));
addpath(genpath('DEA Files')); 
addpath(genpath('Gamma Files')); 

tht0vec=zeros(M,1);
tht0vec(M-nEdge(end)+1:M)=tht0*ones(nEdge(end),1);
tht0vec(1:nEdge(1))=(0.5*pi-tht0)*ones(nEdge(1),1);

stage =1
%% BOUNDARY CONDITIONS AND FFT(BC) 
BC = zeros(M,N);
utildeB = zeros(M,N);
Ftilde_BC=zeros(M,N);
for m=1:NVert
    for mv =1:nEdge(m) 
        mt=mv+cEdge(m);
        for nn=1:N    
            n=nn-1;
            tn = n*dt;
            l=(lambda^n);
%             %BC for wave on left Edge
            if m==NVert || m==1
                BC(mt,nn)=l*(cos(tht0)*nx(mt)+sin(tht0)*ny(mt))*F_BC(xi(mt),yi(mt),tht0,tn,c,alpha,t0); %  d/dn = d/dx
            else
                BC(mt,nn) = 0;
            end 
        end
    end
end
%FFT of Boundary Condition
for m =1:M 
    Ftilde_BC(m,:)=fft(BC(m,:)); 
end

% HYBRID METHOD
%%
stage =2
kdiff = 10000*ones((N/2)+1,1);
utildeDEA=zeros(M,N);
NewTj= zeros(M,1);

%Finds difference between the low frequencies in CQBEM and the switch freq 
Maxkim=max(imag(k));
for ll=1:(N/2)+1
        if imag(k(ll))<0.5*Maxkim
            kdiff(ll) = Freq-abs(real(k(ll)));   
        end
end
%finds the wavenumber value and position in k, for the freq clostest to the switch freq k*
[val,pos] = min(kdiff(kdiff>=0)); 
kfreq=k(pos);

% CALCULATE THE BOUNDARY SOLUTION FOR LOW FREQUENCIES VIA CQBEM
for ll =1:(N/2)+1 
    if abs(real(k(ll))) <=  abs(real(kfreq))  
        utildeB(:,ll) = CQBEM_par_Calc(k(ll),Ftilde_BC(:,ll),M,NVert,nEdge,cEdge,CPi,xi,yi,a,b,xv,yv,CosEdgeAngle,SinEdgeAngle,CL,nx,ny); %CQBEM 
    end
end
utildeB(:,N:-1:(N/2)+2)=conj(utildeB(:,2:(N/2))); %Solutions are complex conjugates 

%%

if DEA==1
    DEAScriptL  %Lshape
else
    %SIMPLE HIGH FREQUENCY APPROXIMATION
    for ll=1:N 
        if abs(real(k(ll))) > abs(real(kfreq))           
            utildeB(:,ll) = Ftilde_BC(:,ll)./(1i*(k(ll))*cos(tht0vec)); %
        end
    end
end
%% INVERSE FFT OF THE BOUNDARY SOLUTION AT ALL FREQUENCY POINTS

Input=zeros(M,N);
uBound=zeros(M,N);
for m=1:M
    Input(m,:) = (ifft(utildeB(m,:))); 
end
%uBound = Boundary solution in time domain.
for n=1:N
    uBound(:,n)=((lambda)^(1-n))*Input(:,n);
end


%% CALCULATE THE INTERIOR SOLUTION
stage =5
Hybrid_Interior_Calc; 

toc

