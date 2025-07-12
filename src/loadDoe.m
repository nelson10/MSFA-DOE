function d = loadDoe(diffractive,N,sigma_d,doe_pitch)
rng("default")
height_max = 1.3745e-6; %2.4e-6;
if(diffractive == 0)
    load('Heightmap.mat');
    d = double(heightmap);
    d = imresize(d,[N N],"box");
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    N1 = 128;
    d = mat2gray(d);
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max; %2.4e-6;
elseif(diffractive == 1)
    % load('MTF_cubic_zeni_-0.3_10.mat'); % 25.5379
    %load('MTF_cubic_zeni_0.1_10_94'); % 25.9567
    %load('MTF_cubic_zeni_0.1_10_best') % 23.8863
    load('MTF_cubic_zeni_0.1_10_previous') % 25.9792
    %load('MTF_cubic_zeni_-0.3_10.mat') % 25.5379
    %load('MTF_cubic_zeni_0.4_40') % 25.4165
    %load('MTF_helical_zeni_0.1_10'); % 27.7628
    %load('MTF_helical_zeni_0.1_10_previous'); % 26.9868
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = DOE2.*A_st;
    N1 = 128;
    d = mat2gray(d);
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max;
elseif(diffractive == 2)
    %Paper: Compact Snapshot Hyperspectral Imaging with Diffracted Rotation
    %Authors: Jeon, Baek, Yi, Dun, Fu, Heidrich, and Kim
    load('spiral_DOE3180_16.mat');
    d = h;
    d = imresize(d,[N N],"nearest");
    d = mat2gray(d).*0.9e-06;
elseif(diffractive == 3)
    alpha = 10*pi; % 5 10 -> inf % 10 d= 1-d -> 0.5 0.8
    M = 1800;
    delta = 3e-6; %pixel size
    [x, y] = meshgrid((-M/2 : M/2-1) *pi*delta);
    d = 1*(x./max(x(:))).^2 + 1*(y./max(y(:))).^2;
    d = mat2gray(mod(alpha.*d, 2*pi));
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = (1-d);
    %imagesc(d)
    %load('fresnel3180_16.mat');
    %d = G3;
    d = imresize(d,[N N],"nearest");
    d = mat2gray(d).*1.3745e-06;%.*0.9e-06;
    d = d.*A_st;
    G3 = d;
    %save('./DOEs/fresnel-DOE_1800.mat','G3')
    %load('Uniform-DOE_1800-fres-0.5-0.8.mat')
elseif(diffractive == 4)
    d = zeros(N,N);
elseif(diffractive == 5)
    %Paper: Extended depth of field through wave-front coding
    % %AuthorS: Edward R. Dowski, Jr., and W. Thomas Cathey
    alpha = 10*pi; %12.5
    M = N;
    delta = doe_pitch; %pixel size
    [x, y] = meshgrid((-N/2 : M/2-1) *pi*delta);
    d = 1*(x./max(x(:))).^3 + 1*(y./max(y(:))).^3;
    d = mat2gray(mod(alpha.*d, 2*pi));
    d = d.*height_max;
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    N1 = 128;
    d = mat2gray(d);
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max;
    %save('./DOEs/Dowski-Cathey_test.mat','d')
    %load('Dowski-Cathey_test.mat')
    %imagesc(d1-mat2gray(d)*height_max)
    %disp('');
elseif(diffractive==6)
    %name = 'Proposed';
    img = imread(fullfile('.', 'DOEs', 'DOE_pattern_159.bmp'));
    if ndims(img) == 3
        img = rgb2gray(img);
    end
    d0 = single(img);
    d0 = imresize(d0, [N, N], 'nearest');
    d = (d0 - min(d0(:))) / (max(d0(:)) - min(d0(:)));
    d = d * height_max;
elseif(diffractive==8)
    load('Phase_delay_profile_6mm.mat');  % no full
    G3 = Phase_delay_profile_6mm_crop_middle;
    G3 = G3(426:end-426,426:end-426);
    d = imresize(G3,[N N],"nearest");
    %load('Wrapped_Phase_delay_profile_6mm_0_2pi.mat'); % no full
    %G3 = Wrapped_Phase_delay_profile_6mm;
    %load('Phase_delay_profile_9.2mm.mat') % full
    %G3 = Phase_delay_profile;
    % load('Wrapped_Phase_delay_profile_9.2mm_0_2pi.mat'); % full
    % G3 = Wrapped_Phase_delay_profile_9mm;
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    N1 = 128;
    d = mat2gray(d);
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max;
elseif(diffractive==7)
    %CGH for generation and focusing of a helical wavefront
    L=1; % number of spirals
    V=0.5;%%Visibility controller
    lambda=0.550*1e-6;%Define wavelength
    f = 11.0;
    del=doe_pitch; % doe pitch
    x=-N/2:N/2-1;
    y=-N/2:N/2-1;
    [X,Y]=meshgrid(x*del,y*del);
    R=sqrt(X.^2+Y.^2);
    alpha = 1;
    betha = 1;
    d = V*exp(alpha*1i*L*(atan2(X,Y))  + betha*1i*((2*pi)/lambda)*(f-sqrt(f*f-R.*R))); %Interference of the object and reference wave
    d = mat2gray(imag(d));
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    N1 = 128;
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max;
    d = d .*A_st;
    G3 = d;
    save("./DOEs/Uniform-DOE_helical_axicon_"+num2str(N),"G3");
end
if(diffractive~=4)
    d = d + normrnd(0,sigma_d,[N,N]);
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
end
