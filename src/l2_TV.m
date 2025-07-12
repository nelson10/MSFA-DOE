function [x_twist] = l2_TV(Y_rgb_noisy,PSF)
x = Y_rgb_noisy;
% remove mean (not necessary)
mu=mean(x(:));
x=x-mu;

B = PSF;
%circularly center
B=fftshift(B);
%normalize
B=B/sum(sum(B));
% convolve
%y = real(ifft2(fft2(B).*fft2(x)));
y = (Y_rgb_noisy);
% set BSNR
BSNR = 15;
Py = var(y(:));
sigma= sqrt((Py/10^(BSNR/10)));

% add noise
%y=y+ sigma*randn(N);

% plot figures
% figure(1); colormap gray;
% imagesc(x); axis off;
% title('Original image')
% figure(2); colormap gray;
% title('Noisy and blurred image')
% imagesc(y); axis off;


% smoothing parameter (empirical setting)
tau = 2e-2*sigma^2/0.56^2;

% extreme eigenvalues (TwIST parameter)
lam1=1e-4;
% TwIST is not very sensitive to this parameter
% rule of thumb: lam1=1e-4 for severyly ill-conditioned% problems
%              : lam1=1e-1 for mildly  ill-conditioned% problems
%              : lam1=1    when A = Unitary matrix


% ------------  TV experiments ---------------------
K=fft2(B);
KC=conj(K);

% handle functions for TwIST
%  convolution operators
A = @(x) real(ifft2(K.*fft2(x)));
AT = @(x) real(ifft2(KC.*fft2(x)));

% denoising function;
tv_iters = 5;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
% TV regularizer;
Phi = @(x) TVnorm(x);
%Phi = @(x) sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));


% start with the wiener filter
varx = var(y(:)); 	% approximate var of x
x0 = real(fft2(KC./(abs(KC).^2+10*sigma^2/varx).*ifft2(y)));

tolA = 1e-4;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function
% falls below 'ToleranceA'
[x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
    TwIST(y,A,tau,...
    'AT', AT, ...
    'lambda',lam1,...
    'True_x', x,...
    'Psi', Psi, ...
    'Phi',Phi, ...
    'Monotone',1,...
    'Initialization',x0,...
    'StopCriterion',1,...
    'ToleranceA',tolA,...
    'Verbose', 1);

end

