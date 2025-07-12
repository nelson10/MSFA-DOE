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
%addpath(genpath('./Kodak24/'));
addpath(genpath('./McMaster/'));
addpath(genpath('./TwIST_v2/'));
addpath(genpath('C:/Users/nelson/OneDrive/Documenten/DeSCI-master-multispectral/dataset'));

a = 1e0; % 1e0 [m]  1e3 [mm]
r = 2.5e-3.*a; % radius of the pupil % 2.5e-3 or 3.0e-3
doe_pitch = 2.0292e-6;%1.86e-6; % DOE pitch
N = round(2*r/doe_pitch);  % Number of grid points per side   %2464
sigma_d = 3e-8;
sigma_s = [0.005, 0.009, 0.015, 0.020];
sigma = sigma_s(2);
algo = 1; %deblurring algorithm, 1 Richarson Lucy, 2 L2-TV
ploton = 1; % To depict groundtruth, mesurement, psf, Doe, Wiener filter recovery
showMTF = 1; % 1 To show MTF otherwise PSF#
dynamicRange = 2^8-1;
alldataset = {'balloons_ms','beads_ms','cd_ms','chart','clay_ms','cloth_ms','egyptian_statue_ms','feathers_ms','flowers_ms','glass_tiles_ms','pompoms_ms','sponges_ms','stuffed_toys_ms','superballs_ms','thread_spools_ms','fake_and_real_beers_ms','face_ms','real_and_fake_peppers_ms','real_and_fake_apples_ms','photo_and_face_ms','paints_ms','oil_painting_ms','jelly_beans_ms','hairs_ms','fake_and_real_tomatoes_ms','fake_and_real_sushi_ms','fake_and_real_strawberries_ms','fake_and_real_peppers_ms','fake_and_real_lemons_ms','fake_and_real_lemon_slices_ms','fake_and_real_food_ms','watercolors_ms'};

mult = 1;
N1 = 512;
L = 31;
[T]=regularMultiplexedSpherePackingCodedAperture(N1,L,mult);
idx = [0 5 6 8 1];
for i=1:5
    diffractive = idx(i); %0 Akpinar TIP-2021, 1 Ours, 2 Spiral-Jeon, 3 Fresnel, 4 without DOE, 5 Dowski-Cathey,6 Oliva,7 Helical axicon, 8 Pinilla et al
    DOE = loadDoe(diffractive,N,sigma_d,doe_pitch);
    %imagesc(DOE),pbaspect([1 1 1]),colormap('jet')
    %imagesc(DOE(426:end-426,426:end-426)),pbaspect([1 1 1])
    %max(DOE(:))
    length(unique(DOE(:)))
    [PSF,deltaS] = computePSF2(DOE);

    p = zeros(32,4);
    s = zeros(32,4);
    sam = zeros(32,4);

    for l=1:32
        dataset = alldataset{l};
        load(dataset);
        I = mat2gray(double(hyperimg))*dynamicRange;
        [Y_md] = DOE_effect(I,PSF,deltaS);
        [Y,Y_orig] = CASSI(Y_md,T,sigma);
        [Xrec] = recover(Y,T);
        %[Xrec] = reconstruction(Y,T,I_orig); % GAP-TV to recover the multispectral datacube
        X = debluring(Xrec,PSF,algo); %Y_orig
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Metrics %%%%%%%%%%%%%%%%%%%%%%%%%%
        disp("Dataset="+num2str(l))
        for k=1:4
            [p(l,k),s(l,k),sam(l,k)] = metrics2(uint8(I),uint8(X(:,:,:,k)));
        end
        %showPlot2(I,X,Y_md,DOE,PSF,l,showMTF,ploton,algo);
    end

    pm = mean(mean(p,2));
    sm = mean(mean(s,2));
    samm = mean(mean(sam,2));

    save("DOE="+num2str(diffractive)+"_Algo="+num2str(algo)+"noise-measure="+num2str(sigma)+'.mat','p','s','sam')
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