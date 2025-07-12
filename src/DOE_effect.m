function [Y_rgb] = DOE_effect(I,PSF,deltaS)
[M,N,L] = size(I);
[Ns1,Ns,L1,K] = size(PSF);
Y_rgb = zeros(Ns1,Ns,L1,K);
X = zeros(Ns1,Ns,L1);
f = zeros(Ns1,Ns,L1);

for l=1:L
    X(:,:,l) =  imresize(I(:,:,l),[Ns1 Ns],'box');
end

for k=1:K
    for l=1:L
        tp = double(X(:,:,l));
        f(:,:,l) = myconv2(tp,PSF(:,:,l,k),deltaS);
    end
    f = norm01(f);
    [Y_rgb(:,:,:,k)] = f;
end
end

function x = norm01(x)
mini=min(x,[],"all");
maxi=max(x,[],"all");
x=(x-mini)/(maxi-mini);
end
