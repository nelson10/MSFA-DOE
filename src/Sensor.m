function [Y,I1,T2,I2] = Sensor(I,T,sigma,shifting)
N1 = size(T,2);
N3 = size(T,1);
K = size(I,3);
D = size(I,4);
dynamicRange = 2^8-1;
rng("default")

for j=1:D
    for i=1:K
        I1(:,:,i,j) = imresize(I(:,:,i,j),[N1 N1],"box");
    end
end

if(shifting==1)
    temp = rand(N3,N3)>=0.5;
    for j=1:D
        for i=1:K
            I2(:,(i-1)+1:(i-1)+N1,i,j) = I1(:,:,i,j);
            T2(:,(i-1)+1:(i-1)+N1,i) = temp;%T(:,:,1);
        end
    end
    [N2,N3,~] =size(T2);
    Y = zeros(N2,N3,D);
else
    Y = zeros(N1,N1,D);
    I2 = I1;
    T2 = T;
end

for j=1:D
    Y1 = I2(:,:,:,j).*T2;
    Y(:,:,j) = mat2gray(sum(Y1,3));
    [m1,n1,~] = size(Y);
    noise = normrnd(0,sigma,[m1,n1]);
    y1 = Y(:,:,j) +  noise;
    y1(y1<=0)=0;
    y1(y1>=1) = 1;
    Y(:,:,j) = y1.*dynamicRange;
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