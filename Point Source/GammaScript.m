
%Find the amplitudes to calculate at k(pos) and k(pos+1)
%VecInf: boundary Density, Tj: No. of directions pointing into the domain along edge j.
VecInf1 = DEA_Calc_PS(dx,k(pos),xv,yv,P0tilde(pos),PSx,PSy,NQ,Omega,MinDamp); %VecInf: boundary Density
[VecInf2,Tj,nSideEls,ThetaLocMat] = DEA_Calc_PS(dx,k(pos+1),xv,yv,P0tilde(pos+1),PSx,PSy,NQ,Omega,MinDamp); %VecInf: boundary Density

 ExtEdge=ones(size(VecInf1));
 % Removes interior edges and takes square root of boundary density to find
 % amplitudes equation (4.8)
 Amp1 = sqrt(VecInf1(ExtEdge>0)); 
 Amp2 = sqrt(VecInf2(ExtEdge>0)); 

Tjp=Tj;
ThetaLocp=squeeze(ThetaLocMat);

CPGlobal = [cEdge,M]; %cummaltive number of elements on edge including total number M
nAmpEdge = nEdge.*Tjp ;  %Number of Amplitudes on each edge 
cAmpEdge = [[0],cumsum(nAmpEdge)]; %Cumlative Number of Amplitudes on each edge 

for i=1:NVert
    for iv=1:nEdge(i)
        it=iv+cEdge(i);
        NewTj(it,1) = Tjp(i);  %number of directions at each collocation point. 
        C(it)= sum(NewTj(1:it)); %Cumalative number of directions at coll pt., used to calculate the global element number in VecInf/AmpTest. 
    end
end
C = [0 C]; %includes the element 0 in the cumalative sum vector

%Creates a table of non-zero amplitudes along each edge / collocation
%point/ direction
ii=zeros(max(Tjp),NVert);
EdgeTab1 = zeros(max(Tj),max(nEdge),NVert);
EdgeTab2 = zeros(max(Tj),max(nEdge),NVert);
for m=1:NVert %Edge index
    for i=1:nEdge(m)
        it=i+cEdge(m); %coll. pt. index
        for n=1:NewTj(it) %Direction index
            if  Amp1(C(it)+n) > max(Amp1)*1e-16
                EdgeTab1(n,i,m) = 1;
            end
            if  Amp2(C(it)+n) > max(Amp2)*1e-16 
                EdgeTab2(n,i,m) = 1;
            end            
        end
    end
%Columns are edges, 0 entries indicate no Amp for that direction at all coll pts along that edge, 
% therefore remove the corresponding gamma coefficent for that direction
    Rowsum1(:,m) = squeeze(sum(EdgeTab1(:,:,m),2));
    Rowsum2(:,m) = squeeze(sum(EdgeTab2(:,:,m),2));

    ii1(:,m) = IdDirRowsum(Rowsum1(1:Tjp(m),m),Tjp(m)); %Entries of ii, indicates the index of directions ThetaLocp that we don't need to solve, for each edge m
    ii2(:,m) = IdDirRowsum(Rowsum2(1:Tjp(m),m),Tjp(m)); %Entries of ii, indicates the index of directions ThetaLocp that we don't need to solve, for each edge m

    Nnzii1(:,m) = nnz(ii1(1:Tjp(m),m)); % number of non zero entries in ii. 
    Nzii1(:,m) = nnz(~ii1(1:Tjp(m),m)); % number of zero entries in ii. 
    Nnzii2(:,m) = nnz(ii2(1:Tjp(m),m)); % number of non zero entries in ii. 
    Nzii2(:,m) = nnz(~ii2(1:Tjp(m),m)); % number of zero entries in ii. 
end

omega1=real(k(pos)); 
omega2=real(k(pos+1));
for m=1:NVert % Loop over edges 
    Indexz = [1:Tjp(m)]'; %Index for gamma for zero ii entries
    Indexnz1 = setdiff(Indexz,ii1(:,m));  %Index for gamma for non zero ii entries
    Indexnz2 = setdiff(Indexz,ii2(:,m));  %Index for gamma for non zero ii entries

    Theta1 = ThetaLocp(m,:);
    Theta2 = ThetaLocp(m,:);
    Theta1(ii1(1:Nnzii1(m),m))=[]; % This removes the directions we are not interested in
    Theta2(ii2(1:Nnzii2(m),m))=[]; % This removes the directions we are not interested in
    
    [B1,Mlhs1] = GammaCalcLin(Amp1(cAmpEdge(m)+1:cAmpEdge(m+1),1),Theta1,Maxutilde1(CPGlobal(m)+1:CPGlobal(m+1)),Tjp(1,m),omega1,CPi(CPGlobal(m)+1:CPGlobal(m+1)),Nzii1(m),nEdge(m),ii1(1:Nnzii1(m),m));
    [B2,Mlhs2] = GammaCalcLin(Amp2(cAmpEdge(m)+1:cAmpEdge(m+1),1),Theta2,Maxutilde2(CPGlobal(m)+1:CPGlobal(m+1)),Tjp(1,m),omega2,CPi(CPGlobal(m)+1:CPGlobal(m+1)),Nzii2(m),nEdge(m),ii2(1:Nnzii2(m),m));
    MLHS{m} = Mlhs2;
    BT{m}=B2;
    GammaVal1 = log(B1)/(1i*omega1);
    GammaVal2 = log(B2)/(1i*omega2);
    
    GammaVecFix1(ii1(1:Nnzii1(m),m),1) = 0; % Reorders Gamma to the correct index position 
    GammaVecFix1(Indexnz1,1) = GammaVal1;

    GammaVecFix2(ii2(1:Nnzii2(m),m),1) = 0; % Reorders Gamma to the correct index position 
    GammaVecFix2(Indexnz2,1) = GammaVal2;

    GammaS1(m,:)=GammaVecFix1.'; 
    GammaS2(m,:)=GammaVecFix2.'; 

end

nGamma=round(real((GammaS2-GammaS1)*omega2*omega1./(2*pi*(omega2-omega1)))); %\nu in thesis.
GammaS=GammaS2+2*pi*nGamma/omega2; %Final Gamma coefficents

% MM=[MLHS{1} zeros(size(Mlhs2)) zeros(size(Mlhs2)) zeros(size(Mlhs2)); zeros(size(Mlhs2)) MLHS{2} zeros(size(Mlhs2)) zeros(size(Mlhs2));zeros(size(Mlhs2)) zeros(size(Mlhs2)) MLHS{3} zeros(size(Mlhs2));zeros(size(Mlhs2)) zeros(size(Mlhs2)) zeros(size(Mlhs2)) MLHS{4}];
% B=[BT{1};BT{2};BT{3};BT{4}];
% 
% PhaseError=norm(MM*B-Maxutilde2)/norm(Maxutilde2)

