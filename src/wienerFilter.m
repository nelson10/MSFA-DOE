%% Created by: Dr. Nelson Eduardo Diaz Diaz
%% Post-doctorado Pontícia Universidad Católica de Valparaíso (PUCV)
%% Date November 11 2023

function [wf] = wienerFilter(PSF,blurred_noisy,delta,noise)
NSR = noise;

H = psf2otf(PSF);
H = fftshift(H);
blurred_noisy = ft2(blurred_noisy,delta);
[~,N,~]=size(blurred_noisy);
H2 = conj(H).*H;
num = conj(H);
den = H2;
wf = ift2(num.*blurred_noisy./(den + NSR),1/(delta)); % 
wf = abs(real(wf));
%wf = mat2gray(wf);
%imagesc(wf);
end