function [vgaptv] = reconstruction(Y,T,orig,shifting)
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

if(shifting==1)
    [M,N1,L,D] = size(vgaptv);
    Xrec = zeros(M,N1-L+1,L,D);
    %ground = zeros(M,N1-L+1,L);
    for l=1:D
        for i=1:L
            Xrec(:,:,i,l) =  vgaptv(:,(i-1)+1:(i-1)+M,i,l);
            %ground(:,:,i) =  orig(:,(i-1)+1:(i-1)+M,i);
            temp = Xrec(:,:,i,l) < 0;
            Xrec(:,:,i,l) = Xrec(:,:,i,l).*~temp;
            Xrec(:,:,i,l) = norm01(Xrec(:,:,i,l));
        end
    end
    vgaptv = Xrec;
end
%save('datos.m',vgaptv);
end

function x = norm01(x)
mini=min(x,[],"all");
maxi=max(x,[],"all");
x=(x-mini)/(maxi-mini);
end