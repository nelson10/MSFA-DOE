function [PSF,h] = createDOI()
%% Parameter setting
d = 1e-3; % diameter 1mm
radii = 0.5e-3;
f = 5e-2; % focal length (50 milimeters)
h00 = 0;%5e-4; % maximum height of the DOE at its center (0.5 milimiters)
L = 25; % Number of spectral channels
NT = 1000;
NR = 1000;
%lbda = round(linspace(lambdaMin,lambdaMax,L));
nabla = 1.52; % refractive index of glass
N = 3; % number of wings
lambdaMin = 4.2* 1e-7; % Minimum visible wavelength
lambdaMax = 6.6* 1e-7; % Maximum visible wavelength
th = linspace(-pi,pi,NT)';
wing = 2*pi/N;
tm = round(NT/N);
x = linspace(-d/2,d/2,NR);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
theta = theta + pi;
Deltah = zeros(NT,NR);
data = zeros(NT*NR,1);
lamb =  zeros(NT,1);
c = 0;
%% Compute radio
for i=1:NR %250
    for j=1:NT %240
        if(r(i,j) <= radii)
            [lb] = lambda(lambdaMin,lambdaMax,N,theta(i,j),wing);
            n = ceil((( (sqrt(r(i,j)^2 + f^2)-f)) )/lb);
            Deltah(i,j) = (n*lb-(sqrt(r(i,j)^2+f^2)-f))./refra(lb*1e6);%nabla;
        end
        c = c + 1;
        data(c,1) = theta(i,j);
        data(c,2) = r(i,j);
        data(c,3) = Deltah(i,j);
    end
end
h = (h00 + Deltah);
PSF = otf2psf(h); % OTF: Optical Transfer Function
PSF = abs(PSF).^2; % PSF: point spread function
%figure,imagesc(h),colormap gray;
%figure,mesh(x,-x,h)


%% Compute lambda(theta)
function [lamb1] = lambda(lambdaMin,lambdaMax,N,theta,wing)
%lambda1 =  zeros(NT,1);
lamb1 = 0;
if(0 <= theta && theta < wing)
    lamb1 = lambdaMin + (lambdaMax - lambdaMin)*(N/(2*pi))*theta;
elseif(theta >= wing)
    %lamb1 = lamb(i-tm);
    lamb1 = lambda(lambdaMin,lambdaMax,N,theta-wing,wing);
    %lamb1 = lambdaMin + (lambdaMax - lambdaMin)*(N/(2*pi))*(theta-wing);
end
end

function[re] = refra(Lb)
IdLens  = 1.5373+0.00829045*Lb^(-2)-0.000211046*Lb^(-4);
re = IdLens - 1;
end
end