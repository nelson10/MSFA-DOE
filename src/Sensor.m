function [Y,T2,I2,detectorName] = Sensor(I,N,sigma,detector)
%% Outputs
% Y -> Measurement
% T2 -> Coded Aperture
% I2 -> Groundtruth

%% Inputs
% I -> Original Imagge
% N -> Spatial resolution
% sigma -> Standard desviation of the noise
% detector -> type of detector

K = size(I,3);
D = size(I,4);
dynamicRange = 2^8-1;
rng("default")

for j=1:D
    for i=1:K
        I1(:,:,i,j) = imresize(I(:,:,i,j),[N N],"box");  % original
    end
end

if(detector==0) % MSFA
    Y = zeros(N,N,D);
    I2 = I1;
    mult = 1;
    %[T]=regularMultiplexedSpherePackingCodedAperture(N,K,mult);
    %save("MSFA.mat","T");
    load("MSFA.mat","T");
    T2 = T;
    detectorName = "MSFA";
elseif(detector==1) % CASSI
    temp = rand(N,N)>=0.5; % Coded aperture with 50% transmittance
    for j=1:D
        for i=1:K
            I2(:,(i-1)+1:(i-1)+N,i,j) = I1(:,:,i,j);
            T2(:,(i-1)+1:(i-1)+N,i) = temp;
        end
    end
    [N2,N3,~] =size(T2);
    Y = zeros(N2,N3,D);
    detectorName = "CASSI";
elseif(detector==2) %Pushbroom
    for j=1:D
        for i=1:K
            noise = normrnd(0,sigma,[N,N]);
            y1 = I1(:,:,i,j) + noise;
            y1(y1<=0)=0;
            y1(y1>=1) = 1;
            I1(:,:,i,j) = y1.*dynamicRange;
        end
    end
    Y = I1;
    T2 = ones(N,N,K);
    I2 = I1;
    detectorName = "Pushbroom";
elseif(detector==3) %RGB
    K = [3 2; 2 1];
    M = kron(K,ones(N/2,N/2));
    for i=1:3
        T(:,:,i) = (M == i);
    end
    %'rggb', 'bggr', 'gbrg', 'grbg'
    Y = zeros(N,N,D);
    I2 = I1;
    T2 = T;
end

if(detector==0 || detector==1 || detector==3) % MSFA, CASSI
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