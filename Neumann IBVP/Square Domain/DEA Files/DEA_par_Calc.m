function [VecInf,Tj,nSideEls,ThetaLocMat]=DEA_par_Calc(dx,waveno,xnodes,ynodes,Force,NQ,t0,Omega)
%Var that can be input.  
%dx = dx; straight swap 
%nSideEls; line 82, 143
%SL; line 144
%stepb; h line 
% maybe just type n=M; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps=1e-15;                          %small parameter near machine precision
epss=1e-8;

%global parameters are abv val of wave vector, damping level and number of vertices of each subdomain

global dAbsP;
global Adamp;
global nVert;  


format long
%parameters
params.omega =real(waveno);% angular freq
params.freq = params.omega/(2*pi);     % frequency                                                
params.rhom = 1;  % density
params.speed=1;
params.dAbsP = (1/params.speed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric data for polygonal domains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nOmega=size(Omega,1);
znodes=zeros(size(xnodes));

XYZ=[xnodes' ynodes' znodes'];
XYZcent=(XYZ(Omega(:,1),:)+XYZ(Omega(:,2),:)+XYZ(Omega(:,3),:)+XYZ(Omega(:,4),:))/4;

%Edge lengths: x,y and z compomnent
Lx=zeros(size(Omega));
Ly=zeros(size(Omega));
Lz=zeros(size(Omega));

%side length and cumulative side length
SL=zeros(size(Omega));
CL=zeros(size(Omega));

% %number of rays pointing into a particular subdomain from a particular edge
Tj=zeros(size(Omega));%sending edge
Ti=zeros(size(Omega));%receiving edge

DiffY=zeros(size(Omega));
MeanX=zeros(size(Omega));

%Define global direction set and elements 
%multiples of 4 give symmetry in each quadrant
N=4*NQ;

 ThetaC=linspace(0,2*pi,N+1);% ray directions
 ThetaAB=0.5*(ThetaC(1:end-1)+ThetaC(2:end)); % ray intervals
 sf=1; %scaling factor

ThetaStep=ThetaAB(2)-ThetaAB(1);
ThetaC=ThetaC(1:N);
ThetaLocMat = zeros(size(Omega,1),size(Omega,2),2*NQ);

% Vals for splitting into boundary elements
n=zeros(nOmega,1);
nSideEls=zeros(size(Omega)); %no. of boundary elements on each edge

stepb=zeros(nOmega,max(nVert));  % boundary element step size on each edge

tval=0.5; %0.5 is val in my derivation: parameter in source and final density later on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:nOmega
    nVert(j)= nnz(Omega(j,:));
    
    %normal vec to each face
    TDNorm(j,:)=cross((XYZ(Omega(j,2),:)-XYZ(Omega(j,1),:)),(XYZ(Omega(j,3),:)-XYZ(Omega(j,1),:))); %normal pointing in z direction
                                                                                                    %not used for 2D probs
    TDNorm(j,:)=TDNorm(j,:)/norm(TDNorm(j,:)); %scales toa vector of length 1.
    

   %vector direction of intersection line between 2 element planes so RotZ
   %and TDNorm coincide: not needed for 2D probs
    UVLp=-[TDNorm(j,2) -TDNorm(j,1) 0];
    if norm(UVLp)>eps
        UVL(j,:)=UVLp/norm(UVLp);
    else
        UVL(j,:)=UVLp;
    end
        
   
    CRA=TDNorm(j,3); %cosine of rotation angle
    SRA=sqrt(1-CRA.^2); %sine of rotation angle
    
    %Defines general rotated x and y axes into plane of element: in 2D only
    %use standard global x and y axes so RotX=(1,0,0), RotY=(0,1,0),
    %RotZ=(0,0,1)
     RotX(j,:)=[CRA+(1-CRA)*UVL(j,1)^2 UVL(j,1)*UVL(j,2)*(1-CRA)+SRA*UVL(j,3) UVL(j,3)*UVL(j,1)*(1-CRA)-SRA*UVL(j,2)];
     RotY(j,:)=[UVL(j,1)*UVL(j,2)*(1-CRA)-SRA*UVL(j,3) CRA+(1-CRA)*UVL(j,2)^2 UVL(j,3)*UVL(j,2)*(1-CRA)+SRA*UVL(j,1)];
     RotZ(j,:)=[UVL(j,1)*UVL(j,3)*(1-CRA)+SRA*UVL(j,2) UVL(j,3)*UVL(j,2)*(1-CRA)-SRA*UVL(j,1) CRA+(1-CRA)*UVL(j,3)^2];
%     
% % %      
    for jj=1:nVert(j)
        jjp=jj+1;
        if jjp>nVert(j)
            jjp=jjp-nVert(j); %vertex numbering modulo nVert. jj =current vertex, jp = next vertex
        end
               
        %Compute the arc lengths, SL is vector of side lengths and CL is vector of cumulative lengths

        Lx(j,jj)=(xnodes(Omega(j,jjp))-xnodes(Omega(j,jj)));               % x-dist btwn vertices
        Ly(j,jj)=(ynodes(Omega(j,jjp))-ynodes(Omega(j,jj)));               % y-dist btwn vertices
        Lz(j,jj)=(znodes(Omega(j,jjp))-znodes(Omega(j,jj)));               % z-dist btwn vertices
        
        SL(j,jj)=sqrt(((Lx(j,jj))^2)+((Ly(j,jj))^2)+((Lz(j,jj))^2));         % Euclidean dist btwn vertices
        CL(j,jj+1)=CL(j,jj)+SL(j,jj);                         % Cumulative lengths
             
        sfactx(j,jj)=Lx(j,jj)/SL(j,jj);
        sfacty(j,jj)=Ly(j,jj)/SL(j,jj); %useful for finding intersections between rays and edges
        sfactz(j,jj)=Lz(j,jj)/SL(j,jj);
    end
    
    dAbsP(j)=params.dAbsP;% %0.5*CL(j,nVert(j)+1) Slowness-vector length = 1/c 
   % dx=CL(j,nVert(j)+1)/Nelt; %elt size
    Adamp(j)=2*imag(waveno); %Adamp = 1; %imaginary part of k=wavenumber
        
    for iv=1:nVert(j)
        nSideEls(j,iv)=max(1,round(SL(j,iv)/dx)); %No. of elements on each edge %nEdge in Jacobs code 
        stepb(j,iv)=SL(j,iv)/nSideEls(j,iv); %step size on each edge   %  h(j)=SL(j)/nEdge(j); in Jacobs Code 
        n(j)=n(j)+nSideEls(j,iv);              %total number of boundary elements    
    end
end

cSideEls=zeros(nOmega,max(nVert)+1); %cum. no. of boundary elements on each edge

for j=1:nOmega
    for iv=1:nVert(j)
        cSideEls(j,iv+1)=cSideEls(j,iv)+nSideEls(j,iv);
    end
end


EdgeOK=(zeros(sum(nVert))); % defines all permissible edge to edge interactions (matrix entry is 1 if possible, 0 if not)
ACheck=(zeros(nOmega,max(nVert))); %checks for colinear edges to stop interactions between them- is 1 unless 2 edges colinear, then is 0
for j=1:nOmega  
    % defines self-subdomain interactions (diagonal blocks)
   EdgeOK(1+sum(nVert(1:j-1)):nVert(j)+sum(nVert(1:j-1)),1+sum(nVert(1:j-1)):nVert(j)+sum(nVert(1:j-1)))=ones(nVert(j))-eye(nVert(j)); 

   for jv=1:nVert(j)     
       jvp=jv+1;     
       while jvp>nVert(j)
             jvp=jvp-nVert(j);
       end
       jvpp=jv+2;
       while jvpp>nVert(j)
             jvpp=jvpp-nVert(j);
       end
       
       PListP(:,jv,j)=[Omega(j,jv); Omega(j,jvp)];
       %changed for mesh case - must change back in 2D general domains
       ACheck(j,jv)=det([1 xnodes(Omega(j,jv)) ynodes(Omega(j,jv));1  xnodes(Omega(j,jvp)) ynodes(Omega(j,jvp));1  xnodes(Omega(j,jvpp)) ynodes(Omega(j,jvpp)) ]);
       %ACheck(j,jv)=1; % simplification for 3D surface mesh 
       
       if ACheck(j,jv)==0 % removes collinear mapping
          
           EdgeOK(jv+sum(nVert(1:j-1)),jvp+sum(nVert(1:j-1)))=0;
           EdgeOK(jvp+sum(nVert(1:j-1)),jv+sum(nVert(1:j-1)))=0;
           
       end
          
   end
end

FreeEdges=ones(nOmega,max(nVert));% edges that don't connnect 2 subdomains

for j=1:nOmega
    for jv=1:nVert(j)
        jvp=jv+1;     
        while jvp>nVert(j)
              jvp=jvp-nVert(j);
        end
        jvm=jv-1;     
        while jvm<1
              jvm=jvm+nVert(j);
        end
        %% Control here to turn curvature reflections off/on
               
        for i=1:nOmega %calculates off-diagonal blocks of EdgeOK
            if j~=i
                for iv=1:nVert(i)
                    if PListP(1:2,jv,j)==PListP(2:-1:1,iv,i) % tests if matching edges for mapping to connected edges                       
                       
                        EdgeOK(iv+sum(nVert(1:i-1)),sum(nVert(1:j-1))+1:sum(nVert(1:j-1))+nVert(j))=ones(1,nVert(j));               
                        EdgeOK(iv+sum(nVert(1:i-1)),jv+sum(nVert(1:j-1)))=0;                    
                        if ACheck(j,jv)==0           %prevents colinear transmission            
                            EdgeOK(iv+sum(nVert(1:i-1)),jvp+sum(nVert(1:j-1)))=0;                             
                        end
                        if ACheck(j,jvm)==0                       
                            EdgeOK(iv+sum(nVert(1:i-1)),jvm+sum(nVert(1:j-1)))=0;                             
                        end                    
                        FreeEdges(j,jv)=0; %sets to zero the non-free edges
                    end
                end
            end
        end
    end
end
%% source Vector data
Source=0; %0 for bc, 1 for point sc
ScDomVert=[];
if Source==0
    % Case 1: BC along edge x=a;
    a=min(xnodes);
    b=min(ynodes);
  %  t0=0*pi;%-0.25*pi;%0;%pi/4;% %source direction - is now an input
%Helmholtz
     %PreFactB=params.rhom/2; % constant prefactor for source term
     
     % Helmholtz Neumann BC
    % waveno=(params.omega/params.speed)+1i*Adamp(1)/2
     
     BCForce = abs(Force).^2;
     PreFactB= BCForce*(1/((cos(t0)*abs(waveno))^2)); % constant prefactor for source term    ???
    
     if t0>0
        PreFactB2= BCForce*(1/((cos(0.5*pi-t0)*abs(waveno))^2)); % constant prefactor for source term    ???
      end

else
    % Case 2: pt source at (x,y)=(Ox,Oy);
    %ScDom=1; % pick a sub-domain for the source point
   
    %final sol coords at all centroids

%     Ox=XYZcent(ScDom,1);    %sets the source point to be the centroid of element nSource... easiest for initial 3d application
%     Oy=XYZcent(ScDom,2);
%     Oz=XYZcent(ScDom,3); 
    
    Ox=XYZ(ScVert,1);    %sets the source point to be the vertex
    Oy=XYZ(ScVert,2);
    Oz=XYZ(ScVert,3); 
    ScDom=[];
    for j=1:nOmega
        for jv=1:nVert(j)
            if ScVert==Omega(j,jv)
               ScDom=[j ScDom];
            end
        end
    end

    ScDom=sort(ScDom);
    for j=1:nOmega
         for jv=1:nVert(j)        
             if any(j==ScDom) && ScVert~=PListP(1,jv,j) && ScVert~=PListP(2,jv,j)
                ScDomVert=[j jv;ScDomVert];
             end
         end
    end
     
    [values, order] = sort(ScDomVert(:,1));
    ScDomVert = ScDomVert(order,:);

    %Ox=0.5;%mean(xnodes(Omega(ScDom,1:nVert(ScDom)))); % put source point in center of chosen subdomain
    %Oy=0.5;%mean(ynodes(Omega(ScDom,1:nVert(ScDom))));
    %Oz=0;    
    PreFact=params.omega*params.rhom/(8*pi);
  
end

%addpath(genpath('distmesh'));
PpOm=0.1; %how fine a triangulation to use for result visualisation

%stage=2
%% Starts main loop over edges (equiv vertices) jv transmitting from 
%Convention: an edge is labelled by its 1st vertex with anti-clockwise orientation
for j=1:nOmega   %loop over subdomains where the ray starts
    for jv=1:nVert(j)    %loop over vertices (equiv edges) in the subdomain j where the ray starts
        
        jvp=jv+1;
         
         while jvp>nVert(j)
            jvp=jvp-nVert(j);
         end
                                   
         xva=xnodes(Omega(j,jv)); %edge start and end points in subdomain j
         xvb=xnodes(Omega(j,jvp));
         yva=ynodes(Omega(j,jv));
         yvb=ynodes(Omega(j,jvp));  
         zva=znodes(Omega(j,jv));
         zvb=znodes(Omega(j,jvp));
         
         EdgeX=dot([xvb-xva yvb-yva zvb-zva],RotX(j,:));
         EdgeY=dot([xvb-xva yvb-yva zvb-zva],RotY(j,:));
        % EdgeZ=dot([xvb-xva yvb-yva zvb-zva],RotZ(j,:));
         
         lowercut=atan2(EdgeY,EdgeX); %minimum angle pointing into Omega
        
         while lowercut<0 % angle range 0 to 2pi
           lowercut=lowercut+2*pi; 
        end
        
        uppercut=lowercut+pi; %maximum angle pointing into Omega%
        
        while uppercut>=2*pi %modulo 2pi
            uppercut=uppercut-2*pi;
        end
        
        if lowercut<uppercut    
            ThetaI=find(ThetaC>=lowercut & ThetaC<uppercut);  %finds index of directions between lowercut and uppercut
            ThetaLoc=0.5*pi-ThetaC(ThetaI)+lowercut; %converts directions ThetaC(ThetaI) to local ones from -pi/2 to pi/2
            
        else     
            ThetaI=[find(ThetaC>=lowercut) find(ThetaC<uppercut)];  %finds index of directions between lowercut and uppercut: 0 to uppercut plus lowercut to 2pi
            ThetaLocL=0.5*pi-ThetaC(ThetaI)+lowercut;
            ThetaLocU=-0.5*pi-ThetaC(ThetaI)+uppercut; %extra complication due to 2 separate ranges 
            ThetaLocL(abs(ThetaLocL)>0.5*pi)=[];
            ThetaLocU(abs(ThetaLocU)>0.5*pi)=[];
            ThetaLoc=[ThetaLocL ThetaLocU]; %converts directions ThetaC(ThetaI) to local ones from -pi/2 to pi/2
           
        end
        
        ThetaLocMat(j,jv,:) = ThetaLoc;
     
        Tj(j,jv)=length(ThetaI); % number of local directions from edge j 
      
        alpha=atan2(EdgeX,-EdgeY); %global direction for normal vector to edge in range -pi to pi.
      
        %block sizes for defining matrix system correspond to number of local directions for each/every element
       
        jBlockA=(sum(sum(Tj(1:j-1,1:end).*nSideEls(1:j-1,1:end)))+sum(Tj(j,1:jv-1).*nSideEls(j,1:jv-1)))+1;
       
        jBlockB=(sum(sum(Tj(1:j-1,1:end).*nSideEls(1:j-1,1:end)))+sum(Tj(j,1:jv).*nSideEls(j,1:jv)));                    
       
        [ajj, bjj, xaj, yaj, zaj, xbj, ybj, zbj] = EltPropsE2(xva,yva,zva,sfactx(j,jv),sfacty(j,jv),sfactz(j,jv),CL(j,jv),stepb(j,jv),nSideEls(j,jv));
        %ajj=arclength values at start of current boundary element that ray
        %is sent from
         %bjj=arclength values at end of current boundary element that ray
         %is sent from
      %xaj,yaj,zaj = cartsian coords on sending element at ajj
        %xbj,ybj,zbj = cartsian coords of sending element at bjj
        
        % for the interior density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PlotCount=0;%zeros(PlotSize(j),1);
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %%%%%%%%%%%SOURCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

        if Source==0
            PreFactBs=zeros(Tj(j,jv),nSideEls(j,jv));         
            % case 1: edge BC along any line with 0 x-coordinate along its length 
            ScVec(jBlockA:jBlockB,1)=zeros(Tj(j,jv)*nSideEls(j,jv),1); % initialise as a vector of zeros
           % midpart=zeros(nSideEls(j,jv),1);
            ScDir=zeros(Tj(j,jv),1);

            if abs(xva-a)<eps && abs(xvb-a)<eps % && abs(ybj)<(0.75+eps) && abs(yaj)>(0.25-eps)
           % if j==1 && jv==4
                % for direction =t0 rads.  
                        

                ScDir(abs(ThetaC(ThetaI)-t0)<=0.51*ThetaStep)=1; % to identify which of global direction set corresponds to the initial density direction t0 
                ScDir(abs(ThetaC(ThetaI)-(t0+2*pi))<=0.51*ThetaStep)=1;% in case t0=0;
                 ScDir=repmat(ScDir,nSideEls(j,jv),1);
              %  midpart(abs(yaj)<(0.75+eps) & abs(ybj)>(0.25-eps))=1; %% nb - yaj>ybj!
              % ScDir=kron(midpart,ScDir);
                  %Made an Edit here, full stop before PreFactB,
                  
                  
                  
                  for loop2=1:nSideEls(j,jv)                  
                      PreFactBs(:,loop2)=ones(Tj(j,jv),1).*PreFactB(loop2+cSideEls(j,jv));
                  end
             
                ScVec(jBlockA:jBlockB,1)=sparse((PreFactBs(:)).*ScDir); %.*(stepb(j,jv)^tval));   %  /SL(j,jv)  
            end
            if t0>0 && abs(yva-b)<eps && abs(yvb-b)<eps
                ScDir(abs(ThetaC(ThetaI)-t0)<=0.51*ThetaStep)=1; % to identify which of global direction set corresponds to the initial density direction t0 
                ScDir(abs(ThetaC(ThetaI)-(t0+2*pi))<=0.51*ThetaStep)=1;% in case t0=0;
                 ScDir=repmat(ScDir,nSideEls(j,jv),1);
              %  midpart(abs(yaj)<(0.75+eps) & abs(ybj)>(0.25-eps))=1; %% nb - yaj>ybj!
              % ScDir=kron(midpart,ScDir);
                  %Made an Edit here, full stop before PreFactB,
                  
                  
                  
                  for loop2=1:nSideEls(j,jv)                  
                      PreFactBs(:,loop2)=ones(Tj(j,jv),1).*PreFactB2(loop2+cSideEls(j,jv));
                  end
             
                ScVec(jBlockA:jBlockB,1)=sparse((PreFactBs(:)).*ScDir); %.*
             end
        else
            %%
        % Case 2: pt source at (x,y)=(Ox,Oy,Oz);
        
            ScInt=zeros(Tj(j,jv),nSideEls(j,jv));
            
            
            % first consider reflected contributions
            % 2nd and 3rd conditions below not needed for source inside domain
            if any(j==ScDom) && ScVert~=PListP(1,jv,j) && ScVert~=PListP(2,jv,j)% restrict to source point domains
            %need to integrate ScTerm entries and place into right positions in ScVec
            
               if FreeEdges(j,jv)==0
                        
                   for loop1=1:Tj(j,jv)
                       for loop2=1:nSideEls(j,jv)                  
                                               
                            SR = @(s) SourceIntegrandR(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(j),ThetaLoc(loop1),ThetaStep);
         
                            ScInt(loop1,loop2)=integral(SR,0,stepb(j,jv));                                     
                                                     
                       end
                   end
               else
                        
                  for loop1=1:Tj(j,jv)
                     for loop2=1:nSideEls(j,jv)                  
                                               
                          SE = @(s) SourceIntegrandE(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(j),ThetaLoc(loop1),ThetaStep);
                                
            
                          ScInt(loop1,loop2)=integral(SE,0,stepb(j,jv));                                     
                     end
                  end
               end
              
            end           
            % now consider transmitted contributions
%             for sv=1:nVert(ScDom) % for source inside domain
              for sv=1:length(ScDom)
               
               %  if PListP(1:2,sv,ScDom)==PListP(2:-1:1,jv,j) %for source in domain:  restrict to source domain connecting edges 
                  if PListP(1:2,ScDomVert(sv,2),ScDom(sv))==PListP(2:-1:1,jv,j)  %for source on vertex
               %                      
               
                    for loop1=1:Tj(j,jv)
                         for loop2=1:nSideEls(j,jv)                  
%                           

                             ST = @(s) SourceIntegrandT(s,xaj(loop2),yaj(loop2),zaj(loop2),xbj(loop2),ybj(loop2),zbj(loop2),Ox,Oy,Oz,Adamp(ScDom(sv)),ThetaLoc(loop1),ThetaStep);
%                          
                             ScInt(loop1,loop2)=integral(ST,0,stepb(j,jv));         
                         end
                    end
                 end
              end

            ScVec(jBlockA:jBlockB,1)=sparse(PreFact.*ScInt(:));
        end
          
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to here
            
           
        for i=1:nOmega   %loop over subdomains where the ray arrives
            for iv=1:nVert(i)    %loop over edges where the ray arrives
                          
                   ivp=iv+1;
         
                    while ivp>nVert(i)
                        ivp=ivp-nVert(i);
                    end
                                
                    xvai=xnodes(Omega(i,iv));
                    xvbi=xnodes(Omega(i,ivp));
                    yvai=ynodes(Omega(i,iv));
                    yvbi=ynodes(Omega(i,ivp));
                    zvai=znodes(Omega(i,iv));
                    zvbi=znodes(Omega(i,ivp));
         
                    TurnAxis=[xvbi-xvai yvbi-yvai zvbi-zvai];
                    
                    EdgeXi=dot(TurnAxis,RotX(i,:));
                    EdgeYi=dot(TurnAxis,RotY(i,:));
                    % EdgeZ=dot([xvb-xva yvb-yva zvb-zva],RotZ(j,:));
         
                    lowercuti=atan2(EdgeYi,EdgeXi);
                            
                    while lowercuti<0
                        lowercuti=lowercuti+2*pi; 
                    end
            
                    uppercuti=lowercuti+pi;
        
                    while uppercuti>=2*pi
                        uppercuti=uppercuti-2*pi;
                    end
        
                    if lowercuti<uppercuti      
                        ThetaIi=find(ThetaC>=lowercuti & ThetaC<uppercuti);
                        ThetaLoci=0.5*pi-ThetaC(ThetaIi)+lowercuti;
            
                    else     
                        ThetaIi=[find(ThetaC>=lowercuti) find(ThetaC<uppercuti)];
                        ThetaLocLi=0.5*pi-ThetaC(ThetaIi)+lowercuti;
                        ThetaLocUi=-0.5*pi-ThetaC(ThetaIi)+uppercuti;
                        ThetaLocLi(abs(ThetaLocLi)>0.5*pi)=[];
                        ThetaLocUi(abs(ThetaLocUi)>0.5*pi)=[];
                        ThetaLoci=[ThetaLocLi ThetaLocUi];
                    end
                    Ti(i,iv)=length(ThetaIi);
                    alphai=atan2(EdgeXi,-EdgeYi);
                    
                    ThetaMat=kron(ones(nSideEls(j,jv),1),ThetaC(ThetaI));
                    ThetaIMat=kron(ones(nSideEls(i,iv),1),ThetaC(ThetaI));
                   
                    
                    iBlockA=(sum(sum(Ti(1:i-1,1:end).*nSideEls(1:i-1,1:end)))+sum(Ti(i,1:iv-1).*nSideEls(i,1:iv-1)))+1;
                    iBlockB=(sum(sum(Ti(1:i-1,1:end).*nSideEls(1:i-1,1:end)))+sum(Ti(i,1:iv).*nSideEls(i,1:iv)));
                    
                 if EdgeOK(iv+sum(nVert(1:i-1)),jv+sum(nVert(1:j-1)))==1 %avoids parallel or non-connecting edges
                    
                   %update
                   RT=ones(Ti(i,iv),Tj(j,jv));
                   %
       
                    StepProd=(stepb(i,iv)^(-0.5))*(stepb(j,jv)^(-0.5));%.*SL(i,iv)/SL(j,jv);
                                                               
                    Bp=zeros(nSideEls(i,iv),Tj(j,jv)*nSideEls(j,jv));

                    [aii, bii, xai, yai, zai, xbi, ybi, zbi] = EltPropsE2(xvai,yvai,zvai,sfactx(i,iv),sfacty(i,iv),sfactz(i,iv),CL(i,iv),stepb(i,iv),nSideEls(i,iv));
          
                      xajMat=(xaj*ones(1,Tj(j,jv)))';
                      yajMat=(yaj*ones(1,Tj(j,jv)))';
                      zajMat=(zaj*ones(1,Tj(j,jv)))';
% % %                    
                       xbjMat=(xbj*ones(1,Tj(j,jv)))';
                       ybjMat=(ybj*ones(1,Tj(j,jv)))';
                       zbjMat=(zbj*ones(1,Tj(j,jv)))';

                       %Another point on ray backprojected from (xai,yai,zai)
                    XRA=xai+cos(ThetaIMat)*RotX(j,1)+sin(ThetaIMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRA=yai+cos(ThetaIMat)*RotX(j,2)+sin(ThetaIMat)*RotY(j,2);
                    ZRA=zai+cos(ThetaIMat)*RotX(j,3)+sin(ThetaIMat)*RotY(j,3);

                    %Another point on ray backprojected from (xbi,ybi,zbi)
                    XRB=xbi+cos(ThetaIMat)*RotX(j,1)+sin(ThetaIMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRB=ybi+cos(ThetaIMat)*RotX(j,2)+sin(ThetaIMat)*RotY(j,2);
                    ZRB=zbi+cos(ThetaIMat)*RotX(j,3)+sin(ThetaIMat)*RotY(j,3);
                   
                    
                    Xa=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Ya=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Za=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Xb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Yb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    Zb=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaAj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbAj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaBj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbBj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DaCj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    DbCj=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                   
                    xcj=0.5*(xaj+xbj);
                    ycj=0.5*(yaj+ybj);
                    zcj=0.5*(zaj+zbj);
                    
                    for ie=1:nSideEls(j,jv)
                        
                        [Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), dummy] = Intersection3D(XRA,YRA,ZRA,xai,yai,zai,xaj(ie),yaj(ie),zaj(ie),xbj(ie),ybj(ie),zbj(ie));
                        [Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv)), dummy] = Intersection3D(XRB,YRB,ZRB,xbi,ybi,zbi,xaj(ie),yaj(ie),zaj(ie),xbj(ie),ybj(ie),zbj(ie));
                   
                        DaAj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xaj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-yaj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zaj(ie)).^2);
                        DbAj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xaj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-yaj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zaj(ie)).^2);
                        DaBj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xbj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ybj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zbj(ie)).^2);
                        DbBj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xbj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ybj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zbj(ie)).^2);
                       
                        DaCj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xa(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xcj(ie)).^2+(Ya(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ycj(ie)).^2+(Za(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zcj(ie)).^2);
                        DbCj(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))=sqrt((Xb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-xcj(ie)).^2+(Yb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-ycj(ie)).^2+(Zb(:,1+(ie-1)*Tj(j,jv):ie*Tj(j,jv))-zcj(ie)).^2);                 
                    
                    end
                    
                    AA=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BB=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BpU=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    BpL=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    
                         XAj=repmat(xajMat(:)',nSideEls(i,iv),1);
                         YAj=repmat(yajMat(:)',nSideEls(i,iv),1);
                         ZAj=repmat(zajMat(:)',nSideEls(i,iv),1);
                         XBj=repmat(xbjMat(:)',nSideEls(i,iv),1);
                         YBj=repmat(ybjMat(:)',nSideEls(i,iv),1);
                         ZBj=repmat(zbjMat(:)',nSideEls(i,iv),1);
                        
                                       
                    AA(DaAj>DaBj)=DaCj(DaAj>DaBj);
                    AA(DaAj<DaBj)=-DaCj(DaAj<DaBj);
                    
                    BB(DbAj>DbBj)=DbCj(DbAj>DbBj);
                    BB(DbAj<DbBj)=-DbCj(DbAj<DbBj);
                    
                  
                    BpU(AA>BB)=bsxfun(@min,0.5*stepb(j,jv),AA(AA>BB));
                    BpL(AA>BB)=bsxfun(@max,-0.5*stepb(j,jv),BB(AA>BB));                       
                    BpU(BB>AA)=bsxfun(@min,0.5*stepb(j,jv),BB(BB>AA));
                    BpL(BB>AA)=bsxfun(@max,-0.5*stepb(j,jv),AA(BB>AA));
           
                    %AA,BB correspond to the back-projection of the ray and 
                    %corresponds to integration limit inside the element,
                    %+ or - 0.5DSj is + or - half element length and
                    %coresponds to (xaj,yaj).
                 
                    if Adamp(j)<eps
                        
                                            
                        Bp(BpU>BpL)=BpU(BpU>BpL)-BpL(BpU>BpL);
                            
                    else
                                               
                        %back proj orientations preserve direction of
                        %destination edge
                        
                         % when preimage outside element, set XBj, YBj to be
                        % lower element bdry and A to be upper
                        
                        
                        if i==j
                             XAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Xb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);
                             YAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Yb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);                           
                             ZAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Zb(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);                           

                             XBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Xa(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);
                             YBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Ya(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);                   
                             ZBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Za(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);                   
                                                                       
                        else
                             XAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Xa(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);
                             YAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Ya(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);                           
                             ZAj(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps)=Za(BpL<BpU & abs(BpL+0.5*stepb(j,jv))>eps);                           

                             XBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Xb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);
                             YBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Yb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);                   
                             ZBj(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps)=Zb(BpL<BpU & abs(BpU-0.5*stepb(j,jv))>eps);                   
                        end
                   GlobMat=repmat(ThetaIMat,1,nSideEls(j,jv));
 
                        %Another point on ray backprojected from (xai,yai,zai)
                    XRAi=XAj+cos(GlobMat)*RotX(j,1)+sin(GlobMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRAi=YAj+cos(GlobMat)*RotX(j,2)+sin(GlobMat)*RotY(j,2);
                    ZRAi=ZAj+cos(GlobMat)*RotX(j,3)+sin(GlobMat)*RotY(j,3);
                    
                    %Another point on ray backprojected from (xbi,ybi,zbi)
                    XRBi=XBj+cos(GlobMat)*RotX(j,1)+sin(GlobMat)*RotY(j,1); %Another point on the ray moving along a unit vector in direction DirEdge
                    YRBi=YBj+cos(GlobMat)*RotX(j,2)+sin(GlobMat)*RotY(j,2);
                    ZRBi=ZBj+cos(GlobMat)*RotX(j,3)+sin(GlobMat)*RotY(j,3);
                    
                    XAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    YAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    ZAi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    XBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    YBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));
                    ZBi=zeros(nSideEls(i,iv),nSideEls(j,jv)*Tj(j,jv));

                    for ie=1:nSideEls(i,iv)
                        [XAi(ie,:), YAi(ie,:), ZAi(ie,:), dummy] = Intersection3D(XRAi(ie,:),YRAi(ie,:),ZRAi(ie,:),XAj(ie,:),YAj(ie,:),ZAj(ie,:),xai(ie),yai(ie),zai(ie),xbi(ie),ybi(ie),zbi(ie));
                        [XBi(ie,:), YBi(ie,:), ZBi(ie,:), dummy] = Intersection3D(XRBi(ie,:),YRBi(ie,:),ZRBi(ie,:),XBj(ie,:),YBj(ie,:),ZBj(ie,:),xai(ie),yai(ie),zai(ie),xbi(ie),ybi(ie),zbi(ie));
                    end
                        % for paralel edges will have problem with La=Lb
                       
                        
                        Lb=sqrt((XBi-XBj).^2+(YBi-YBj).^2+(ZBi-ZBj).^2);
                        La=sqrt((XAi-XAj).^2+(YAi-YAj).^2+(ZAi-ZAj).^2);
                        dLds=(Lb-La)./(BpU-BpL);
                              
                        if La(BpU>BpL)~=Lb(BpU>BpL)
                                                        
                            Bp(BpU>BpL)=(-exp(-Adamp(j)*Lb(BpU>BpL))+exp(-Adamp(j)*La(BpU>BpL)))./(dLds(BpU>BpL)*Adamp(j));                       
                                                        
                        else %parallel edge case, L is const
                            
                            Bp(BpU>BpL)=(exp(-Adamp(j)*Lb(BpU>BpL))).*(BpU(BpU>BpL)-BpL(BpU>BpL));

                        end
                       
                              
                    end        
                    
                  
                    ThetaCMat=kron(ThetaC(ThetaIi)',ones(1,Tj(j,jv)));
                                           
                    DirMap=zeros(Ti(i,iv),Tj(j,jv));
                   
                     if i==j

                        AngleRefIn=pi+ThetaC(ThetaI);
                        AngleRefOut=mod(2*alphai-AngleRefIn,2*pi);          
                        AngleRefM=kron(ones(Ti(i,iv),1),AngleRefOut);
                                          
                        DirMap(abs(AngleRefM-ThetaCMat)<0.5*ThetaStep)=1;%
                        
                     
                        DirMap(mod(abs(ThetaCMat-alphai+0.5*pi),2*pi)==0)=0;
                        
%                         ThetaIIMat=kron(ones(nSideEls(i,iv),1),ThetaC(ThetaI));
                        ThetaIIMat=kron(ones(Ti(i,iv),1),ThetaC(ThetaI));

                        DirMap(mod(abs(ThetaIIMat-alpha+0.5*pi),2*pi)==0)=0;
                      
                        if FreeEdges(i,iv)==0           %reflection law                       
                           RT=0*ones(Ti(i,iv),Tj(j,jv));
                        else             %wall damping factor          
                           RT= 1*ones(Ti(i,iv),Tj(j,jv));  
                        end                         
                      
                    else
                       
                        TurnAng=acos(dot(TDNorm(i,:),TDNorm(j,:)));
                                  
                        TurnAxisU=-TurnAxis./norm(TurnAxis);
                     
                        ConTest=dot(XYZcent(i,:)-XYZcent(j,:),TDNorm(i,:)-TDNorm(j,:));
                        %Tests if faces convex of concave so rotate about
                        %correct angle
                        
                        if ConTest<0 % correction for concave edges
                            TurnAng=-TurnAng;
                        end
                        
                        cT=cos(TurnAng);
                        sT=sin(TurnAng);
                        %TurnX is the local axis in element j, rotated into the plane of element i   
                        TurnX=[cT+(1-cT)*TurnAxisU(1)^2 (1-cT)*TurnAxisU(1)*TurnAxisU(2)-sT*TurnAxisU(3) (1-cT)*TurnAxisU(1)*TurnAxisU(3)+sT*TurnAxisU(2); (1-cT)*TurnAxisU(1)*TurnAxisU(2)+sT*TurnAxisU(3) cT+(1-cT)*TurnAxisU(2)^2 (1-cT)*TurnAxisU(2)*TurnAxisU(3)-sT*TurnAxisU(1); (1-cT)*TurnAxisU(1)*TurnAxisU(3)-sT*TurnAxisU(2) (1-cT)*TurnAxisU(2)*TurnAxisU(3)+sT*TurnAxisU(1) cT+(1-cT)*TurnAxisU(3)^2]*(RotX(j,:)');
                        
                        SignTest=dot(TurnX,RotY(i,:));
                       
                        % Rotates the directions in element j to orient with the
                        % local x-axis in element i
                        AngleTran=mod(ThetaC(ThetaI)+sign(SignTest)*abs(acos(min(dot(RotX(i,:),TurnX),1))),2*pi);%
                         
                        AngleTran(abs(AngleTran-2*pi)<2e-8)=0;%AngleTran(abs(AngleTran-2*pi)<2e-8)-AngleTran(abs(AngleTran-2*pi)<2e-8);
                       % end
                        AngleTranM=kron(ones(Ti(i,iv),1),AngleTran);
                                            
                        DirMap(abs(AngleTranM-ThetaCMat)<0.5*ThetaStep)=1;  %% needs the be corrected fo matrices always same size
                        
                        DirMap(abs(AngleTranM-ThetaCMat-2*pi)<0.5*ThetaStep)=1; %ADDED SINCE NEEDS TO BE MOD 2PI!! %% needs the be corrected fo matrices always same size
                     end
                    
                     Bp(isnan(Bp))=0;
                     Bp=repelem(Bp,Ti(i,iv),1);                    
                     RTE=repmat(RT,nSideEls(i,iv),nSideEls(j,jv));

                     DirMapE=repmat(DirMap,nSideEls(i,iv),nSideEls(j,jv));
                     Bp=RTE.*Bp.*DirMapE;
                    %   
                    TrOpLoc(iBlockA:iBlockB,jBlockA:jBlockB)=StepProd*sparse(Bp);
                   
               end
                

            end
        end
        
    end
    
end
   %%
%stage=3
NM=length(TrOpLoc);
LHS=speye(NM)-TrOpLoc;  % calculate I-T = LHS
VecInf=LHS\(sf*ScVec);  % evaluate LHS^{-1} times rho_0
% plot((abs(ScVec)))
