function [T2]=regularMultiplexedSpherePackingCodedAperture(N,NF,mult)
addpath(genpath('./utils'));
G = zeros(NF,NF,mult);
T = zeros(NF,NF,NF);
T1 = zeros(NF,NF,NF,'logical');
T2 = zeros(N,N,NF);
density = zeros(1,NF-1);
diameter = zeros(1,NF-1);
cn= zeros(1,NF-1);

randomNumber  = 8;
%% Load SOTA Random Coded Aperture
rng(randomNumber);
[~,~,~,G2]=DDDRSNNPLattice(NF,NF);
G(:,:,1) = G2;
Dens_best = 0;
for j=1:500
    T1(:,:,1) = false;
    idx = randperm(NF,mult);
    for i=1:mult
        T1(:,:,1) =  T1(:,:,1)  | G(:,:,1)==idx(i);
    end
    [density1,~]= ComputeDensityIrregularSP(1.*(T1(:,:,1)));
    if(Dens_best < density1)
        best_idx = idx;
        Dens_best = density1;
    end
end

idx = best_idx;
T1(:,:,1) = false;
for i=1:mult
    T1(:,:,1) =  T1(:,:,1)  | G(:,:,1)==idx(i);
end
imagesc(T1(:,:,1))

T(:,:,1) = T1(:,:,1);

for j=1:NF-1
    T(:,:,:)=0;
    T(:,:,1) = T1(:,:,1);
    for i=1:NF-1
        T(:,:,i+1) = circshift(T(:,:,i)',[0 j])';
    end
    unif = sum(T,3);
    if(std(unif(:))==0 && unique(unif(:))==mult)
        [density(j),diameter(j)]= ComputeDensityIrregularSP(1.*T);
        [cn(j)] = ComputeConditionNumber(1.*(T));
    else
        density(j)=0;
        cn(j) = 0;
    end
end
[~,id]=max(density);

for i=2:NF
    T(:,:,i) = circshift(T(:,:,i-1)',id)';
end

r = ceil(N/NF);
B = ones(r,r);

%msfaSize = length(G)*length(B);
%T2 = zeros(msfaSize,msfaSize,mult);
for i=1:NF
    tempo(:,:,i) = kron(B,T(:,:,i));
end
T2 = tempo(1:N,1:N,:);
end