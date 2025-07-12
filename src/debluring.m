function [X_unblur] = debluring(X,PSF,algo)
[~,Ns,L,K] = size(X);
PSF1 = zeros(Ns,Ns,L,K);
X_unblur = zeros(Ns,Ns,L,K);
dynamicRange = 2^8-1;

for k=1:K
    for l=1:L
        PSF1(:,:,l,k) = imresize(PSF(:,:,l,k),[Ns Ns],"box");
    end
end

for k=1:K
    for l=1:L
        y1 = X(:,:,l,k);
        %y1 = mat2gray(y1;
        if (algo==1)
            %if(sigma == 0)
            %   [X_unblur(:,:,l,k)] = wienerFilter(PSF1(:,:,l,k),y1,deltaS,estimated_nsr);
            %else
            %X_unblur(:,:,l,k) = deconvwnr(y1,PSF1(:,:,l,k),R1);
            %X_unblur(:,:,l,k) = deconvreg(y1,PSF1(:,:,l,k),sigma);
            X_unblur(:,:,l,k) = deconvlucy(y1,PSF1(:,:,l,k),15);
            %end
        elseif(algo==2)
            [X_unblur(:,:,l,k)] = (l2_TV(y1,PSF1(:,:,l,k)));
        end
    end
end
if (algo==1)
    X_unblur = mat2gray(X_unblur);
elseif(algo==2)
    X_unblur(X_unblur<=0)=0;
    X_unblur(X_unblur>=1)=1;
end
X_unblur = uint16(X_unblur.*dynamicRange);
end