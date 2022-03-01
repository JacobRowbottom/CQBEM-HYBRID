%Find the amplitudes to calculate at k(pos) and k(pos+1)
%VecInf: boundary Density, Tj: No. of directions pointing into the domain along edge j.

VecInf1= DEA_par_Calc_L(dx,k(pos),xv,yv,Ftilde_BC(:,pos),NQ,tht0,Omega,M,cEdge); %VecInf: boundary Density, Tj: No. of directions pointing into the domain along edge j.
[VecInf2,Tj,nSideEls,ThetaLocMat]= DEA_par_Calc_L(dx,k(pos+1),xv,yv,Ftilde_BC(:,pos+1),NQ,tht0,Omega,M,cEdge); %VecInf: boundary Density, Tj: No. of directions pointing into the d

% For DEA multi-domains (nOmega>1) need to remove entries in VecInf corresponding to "internal" boundaries. 

% For L-shape internal edges are edge 2/4 in subdomain 1 and edge 4/4 in subdomain 2. 
%Create a vector of 1s and 0s with 0 for internal edge coeffs 
 ExtEdge=ones(size(VecInf1));
 ExtEdge(Tj(1,1)*nSideEls(1,1)+1:sum(Tj(1,1:2).*nSideEls(1,1:2)))=zeros(size(Tj(1,2).*nSideEls(1,2)));
 ExtEdge(sum(Tj(1,1:end).*nSideEls(1,1:end))+sum(Tj(2,1:3).*nSideEls(2,1:3))+1:sum(Tj(1,1:end).*nSideEls(1,1:end))+sum(Tj(2,1:end).*nSideEls(2,1:end)))=zeros(size(Tj(2,4).*nSideEls(2,4)));

 Amp1 = sqrt(VecInf1(ExtEdge>0)); % Removes interior edges
 Amp2 = sqrt(VecInf2(ExtEdge>0)); % Removes interior edges

 Tjp=[Tj(1,1); Tj(1,3:4)'; Tj(2,1:3)'];
 ThetaLocp=[ThetaLocMat(1,1,:); permute(ThetaLocMat(1,3:4,:),[2,1,3]); permute(ThetaLocMat(2,1:3,:),[2,1,3])];
 
% Reorder edges from DEA 1,5,6,2,3,4 to BEM 1,2,3,4,5,6
 Amp1Order1 =Amp1(1:Tj(1,1)*nSideEls(1,1)); %edge 1
 Amp1Order3 =Amp1(Tj(1,1)*nSideEls(1,1)+1:sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2));
 Amp1Order2 =Amp1(sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2)+1:end);
% 
 Amp2Order1 =Amp2(1:Tj(1,1)*nSideEls(1,1)); 
 Amp2Order3 =Amp2(Tj(1,1)*nSideEls(1,1)+1:sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2));
 Amp2Order2 =Amp2(sum(Tj(1,1:end).*nSideEls(1,1:end))-Tj(1,2)*nSideEls(1,2)+1:end);

 Amp1Final=[Amp1Order1;Amp1Order2;Amp1Order3]; % Amplitudes at freq switch used to determine gamma
 Amp2Final=[Amp2Order1;Amp2Order2;Amp2Order3]; % Amplitudes at freq switch used to determine gamma

CPGlobal = [cEdge,M]; %cummaltive number of elements on edge including total number M
nAmpEdge = nEdge.*Tjp' ;  %Number of Amplitudes on each edge 
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
EdgeTab1 = zeros(max(max(Tj)),max(nEdge),NVert);
EdgeTab2 = zeros(max(max(Tj)),max(nEdge),NVert);
for m=1:NVert    %Edge index
%     EdgeTab = zeros(Tj(m),nEdge(m),m);% Table which has a 1 in entry if Coll pt. has a non zero amp
    for i=1:nEdge(m)%coll. pt. index
        it=i+cEdge(m);
        for n=1:NewTj(it) %Direction index
            if  Amp1Final(C(it)+n) > max(Amp1Final)*1e-16 
                EdgeTab1(n,i,m) = 1;
            end
            if  Amp2Final(C(it)+n) > max(Amp2Final)*1e-16 
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
    
    [B1,Mlhs1] = GammaCalcLin(Amp1Final(cAmpEdge(m)+1:cAmpEdge(m+1),1),Theta1,Maxutilde1(CPGlobal(m)+1:CPGlobal(m+1)),Tjp(m,1),omega1,CPi(CPGlobal(m)+1:CPGlobal(m+1)),Nzii1(m),nEdge(m),ii1(1:Nnzii1(m),m));
    [B2,Mlhs2] = GammaCalcLin(Amp2Final(cAmpEdge(m)+1:cAmpEdge(m+1),1),Theta2,Maxutilde2(CPGlobal(m)+1:CPGlobal(m+1)),Tjp(m,1),omega2,CPi(CPGlobal(m)+1:CPGlobal(m+1)),Nzii2(m),nEdge(m),ii2(1:Nnzii2(m),m));
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