function [PSF,deltaS] = computePSF(DOE)
a = 1e0; % 1e0 [m]  1e3 [mm]
%z = [0.5 0.8 1.2 inf].*a; %
z = [0.5 1.0 2.0 inf].*a; %
sigma_s = [0.005, 0.009, 0.015, 0.020];
id_sigma = 1;
gamma = 0.125;
algo = 2; % 1 Wiener, 2 l2-TV
ploton = 1; % To depict groundtruth, mesurement, psf, Doe, Wiener filter recovery
showMTF = 1; %0 to show PSF, 1 to show MTF
r = 2.5e-3.*a; % radius of the pupil % 2.5e-3 or 3.0e-3
doe_pitch = 2.0292e-6;%1.86e-6; % DOE pitch
N = round(2*r/doe_pitch);  % Number of grid points per side   %2464
Ns = round(N/2);
M = N;
rng(0);

wln = [6.11e-7 5.43e-7 4.82e-7].*a;  %wavelength [m]
central_wavelen = 5.43e-7;
%zfc = [1.62 1.3 1.04];
L = size(wln,2); % Number of channels
n_lambda = [1.457 1.460 1.463];  % wavelength-dependent refractive index
%zi = 0.0369.*a;  % focal lenght from CMOS to DOE
zi = 2.*r.*a; % The sensor-to-lens distance

%z = [0.025 0.042 0.050 0.065].*a;
%f_lambda=[0.036078, 0.035882, 0.035636];% 36.078, 35.882, 35.636  wavelength-dependent focal lenght
R = 2.3e-3;% %2.3 mm Radius of curvature from lens f_lambda = 5 mm
f_lambda = R ./(n_lambda-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = size(z,2); % Number of depths
h1 = zeros(Ns,Ns,L);
p = zeros(K,1);
si = zeros(K,1);
pmean = zeros(24,1);
smean = zeros(24,1);
X = zeros(Ns,Ns,L,'uint8');
f = zeros(Ns,Ns,L);
Y_rgb = zeros(Ns,Ns,L,K);
wftempo = zeros(Ns,Ns,L);
wf = zeros(Ns,Ns,K,L);
color={'r','g','b'};
H = zeros(Ns,Ns,L);
H1 = zeros(Ns,Ns,L);
Psi = zeros(L,1);
temp = wln*1e9;
wavelength = string(temp);
Y_rgb_noisy = zeros(Ns,Ns,L,K);
h_lambda_z = zeros(N,M,L,K);
I_s = zeros(N,M,L,K);
algorithm = ["Wiener","l2-TV"];
% indices of grid points in aperture
im = ones(M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSF computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:K
    for l=1:L
        Psi(l,k) = (pi/wln(l)).*((1/z(k)) +(1/zi) - (1/f_lambda(l))) .* r^2; % defocus coefficent
    end
end
MaxDefocus = floor(max(abs(Psi(:))));
%% Pupil
deltaS =  1/7*(r*pi)/(8*MaxDefocus); %13 -> best %Eq. 10 TIP-2021  MaxDefocus =21e-6 micro -> zi = 0.0369

[x, y] = meshgrid((-N/2 : M/2-1) .* doe_pitch); %1e-6
A_st = elliptical_crop(im,1)>0;
S = (x.^2 + y.^2)/(r^2);
%% Load Diffractive Optical Element
d = DOE;
for k=1:K
    for l=1:L
        k1 = 2*pi/wln(l); % wavenumber
        Phi = k1*(n_lambda(l)-1).*double(d); %diffractive
        Psi(l,k) = (pi/wln(l)).*((1/z(k)) +(1/zi) - (1/f_lambda(l))) .* r^2; % defocus coefficent
        Q_lamda_z = A_st .* exp(1i*Phi + 1i*Psi(l,k)*S);  %generalized pupil function
        h_lambda_z(:,:,l,k) = abs(ft2(Q_lamda_z, deltaS)).^2; % depth-dependent PSF
        %I_s(:,:,c,k) = mat2gray(myconv2(I(:,:,c),h_lambda_z(:,:,c),delta));
    end
end
%PSF = h_lambda_z;
PSF = imresize(h_lambda_z,[Ns Ns],"box");
end