%% Pontificia Universidad Católica de Valparaíso
%% Author: Dr. Nelson Eduardo Diaz Diaz
%% Simulator for incoherent extended-depth-of-field
%% Date: July 7, 2025
%% Valparaíso, Chile
%% Add path to the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Path definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;

addpath(genpath('./DOEs'));
addpath(genpath('./src'));
addpath(genpath('./reconst'));
addpath(genpath('./GAP-TV'));
addpath(genpath('./utils2'));
addpath(genpath('./RGB'));
addpath(genpath('./metrics'));
addpath(genpath('./Kodak24/'));
addpath(genpath('./McMaster/'));
addpath(genpath('./TwIST_v2/'));
addpath(genpath('./dataset/'));

a = 1e0; % 1e0 [m]  1e3 [mm]
r = 2.5e-3.*a; % radius of the pupil % 2.5e-3 or 3.0e-3
doe_pitch = 2.0292e-6;%1.86e-6; % DOE pitch
N = round(2*r/doe_pitch);  % Number of grid points per side   %2464
sigma_d = 3e-8;
sigma_s = [0 0.005, 0.009, 0.015, 0.020];
sigma = sigma_s(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detector = 2; % 0 MSFA, 1 CASSI, 2 Pushbroom, 3 RGB
algo = 1; %deblurring algorithm, 1 Richarson Lucy, 2 L2-TV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ploton = 0; % To depict groundtruth, mesurement, psf, Doe, Wiener filter recovery
showMTF = 1; % 1 To show MTF otherwise PSF#
dynamicRange = 2^8-1;

N1 = 512;
L = 31;
idx = [1 0 5 8 6 4];%[7 1 0 5 6 8];

for i=1:6
    diffractive = idx(i); %0 Akpinar TIP-2021, 1 Diaz, 2 Spiral-Jeon, 3 Fresnel, 4 without DOE, 5 Dowski-Cathey,6 Oliva,7 Helical axicon, 8 Pinilla et al
    [DOE,doeName] = loadDoe(diffractive,N,sigma_d,doe_pitch);
    %imagesc(DOE),pbaspect([1 1 1]),colormap('jet')
    %imagesc(DOE(426:end-426,426:end-426)),pbaspect([1 1 1])
    %max(DOE(:))
    length(unique(DOE(:)))
    [PSF,deltaS] = computePSF2(DOE,L);
    Ns = size(PSF,2);

    K = size(PSF,4);
    p = zeros(32,K);
    s = zeros(32,K);
    sam = zeros(32,K);

    for l=1:32
        [I] = loadImage(detector,l,dynamicRange);
        [Y_md] = DOE_effect(I,PSF,deltaS);
        [Y,T1,I2,detectorName] = Sensor(Y_md,N1,sigma,detector);  %0 MSFA, 1 CASSI, 2 Pushbroom, 3 RGB
        [Xrec] = reconstruction(Y,T1,I2,detector);
        X = debluring(Xrec,PSF,algo); %Y_orig
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics %%%%%%%%%%%%%%%%%%%%%%%%%%
        disp("Dataset="+num2str(l))
        I = imresize(I,[Ns Ns],"box");
        K = size(p,2);
        for k=1:K
            [p(l,k),s(l,k),sam(l,k)] = metrics2(uint8(I),uint8(X(:,:,:,k)));
        end
        %showPlot2(I,X,Y_md,DOE,PSF,l,showMTF,ploton,algo);
    end

    pm = mean(mean(p,2));
    sm = mean(mean(s,2));
    samm = mean(mean(sam,2));

    save("DOE=_"+num2str(diffractive)+doeName+"_AlgoDebluring=_"+num2str(algo)+"_Detector=_"+detectorName+"_noise-measure=_"+num2str(sigma)+'.mat','p','s','sam')
    dist = ["z1", "z2", "z3", "z4", "z5", "z6", "avg"];
    PSNR_mean = mean(p, 1);
    PSNR_std  =  std(p, 1);
    SSIM_mean = mean(s, 1);
    SSIM_std  =  std(s, 1);
    SAM_mean = mean(sam, 1);
    SAM_std  =  std(sam, 1);
    avg_PSNR = mean(p, 'all');
    avg_SSIM = mean(s, 'all');
    avg_SAM  = mean(sam, 'all');
    std_PSNR = mean(PSNR_std);
    std_SSIM = mean(SSIM_std);
    std_SAM  = mean(SAM_std);
    varNames = {'dist', 'PSNR_avg', 'PSNR_std', 'SSIM_avg', 'SSIM_std', 'SAM_avg', 'SAM_std'};
    Tabla = table(dist', [PSNR_mean avg_PSNR]',[PSNR_std std_PSNR]' , [SSIM_mean avg_SSIM]', [SSIM_std std_SSIM]', [SAM_mean avg_SAM]', [SAM_std std_SAM]', 'VariableNames',varNames )
    writetable(Tabla, names_sensor(shifting+1) +"_"+names_DOE(DOE_idx)+"_"+Algorithm(algo)+"_sigmaID_"+num2str(sigma_id)+"_full.csv")
end

% X = double(X);
%
% for i=1:4
%     RGB = RGB_test(Y_md(:,:,:,i));
%     subplot(2,4,i),imagesc(RGB),title("Blur scr "+num2str(i));
%     RGB = RGB_test(X(:,:,:,i));
%     meanp = mean(p(:,i));
%     subplot(2,4,4+i),imagesc(RGB),title("UnBlur scr "+num2str(i)+"PSNR "+num2str(meanp));
% end