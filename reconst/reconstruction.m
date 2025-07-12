function [vgaptv] = reconstruction(Y,T,orig)
D = size(Y,3);

para.nframe =   1; % number of coded frames in this test
para.MAXB   = 255;
nframe = para.nframe; % number of coded frames in this test
MAXB = para.MAXB;
% [2] apply GAP-Denoise for reconstruction
para.Mfunc  = @(z) A_xy(z,T);
para.Mtfunc = @(z) At_xy_nonorm(z,T);

para.Phisum = sum(T.^2,3);
para.Phisum(para.Phisum==0) = 1;
% common parameters
para.lambda   =     1; % correction coefficiency
para.acc      =     1; % enable GAP-acceleration
para.flag_iqa = false; % disable image quality assessments in iterations

%% [2.1] GAP-TV, ICIP'16
para.denoiser = 'tv'; % TV denoising
para.maxiter  = 300; % maximum iteration
para.tvweight =  15; % weight for TV denoising
para.tviter   =  10; % number of iteration for TV denoising

parfor l=1:D
    [vgaptv(:,:,:,l),psnr_gaptv,ssim_gaptv,tgaptv] = gapdenoise_cacti(T,Y(:,:,l),orig(:,:,:,l),[],para);
end
%save('datos.m',vgaptv);
end