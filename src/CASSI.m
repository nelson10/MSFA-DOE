function [Y,I1] = CASSI(I,T,sigma)
N1 = size(T,2);
K = size(I,3);
D = size(I,4);
Y = zeros(N1,N1,D);
dynamicRange = 2^8-1;
rng("default")

for j=1:D
    for i=1:K
        I1(:,:,i,j) = imresize(I(:,:,i,j),[N1 N1],"box");
    end
end

for j=1:D
    Y1 = I1(:,:,:,j).*T;
    Y(:,:,j) = (sum(Y1,3)).*dynamicRange;
    noise = normrnd(0,sigma,[N1,N1]);
    y1 = Y(:,:,j) +  noise;
    y1(y1<=0)=0;
    y1(y1>=dynamicRange) = dynamicRange;
    Y(:,:,j) = y1;
end

% noise = normrnd(0,sigma,[N1,N1]);
% %noise1 = awgn(ones(N,N),snr,'measured')-1; % noise
% Sn = abs(fft2(noise)).^2; % noise power spectrum
% nA = sum(Sn(:))/numel(noise); % noise average power
% Sf = abs(fft2(f)).^2; % image power spectrum
% fA = sum(Sf(:))/numel(f); % image average power
% R1 = abs(nA) / abs(fA);
% y1 = f +  noise;
%estimated_nsr = sigma / var(y1(:));

end