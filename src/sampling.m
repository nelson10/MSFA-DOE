function [pm,sm] = sampling(I,PSF,deltaS)
L = size(I,3);
N = size(PSF,2)*2;
K = size(PSF,4);
Ns = round(N/2);
ploton = 0;
sigma_s = [0.005, 0.009, 0.015, 0.020];
id_sigma = 1;
algo = 1; % 1 Wiener, 2 l2-TV

for l=1:L
    if(N == 500)
        X(:,:,l) =  I(:,:,l);
    else
        X(:,:,l) =  imresize(I(:,:,l),[Ns Ns]);
    end
end
RGB_X = X;
for k=1:K
    for l=1:L
        tp = double(X(:,:,l));
        f(:,:,l) = myconv2(tp,PSF(:,:,l,k),deltaS);
    end
    f = norm01(f).*255;
    [Y_rgb(:,:,:,k)] = f;%SimulateRGB(f,L1,WlTop,Wlbtn);
end
if(ploton==1)
    figure(1);
    subplot(3,K+1,1),imagesc(X),title("Groundtruth Image "+num2str(j)),pbaspect([1 1 1]),axis off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wiener Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_noise = sigma_s(id_sigma);
for k=1:K
    Y_mesu = Y_rgb(:,:,:,k);
    parfor l=1:L
        %f = Y_mesu(:,:,l);
        noise = normrnd(0,sigma_noise,[Ns,Ns]);
        %noise1 = awgn(ones(N,N),snr,'measured')-1; % noise
        Sn = abs(fft2(noise)).^2; % noise power spectrum
        nA = sum(Sn(:))/numel(noise); % noise average power
        Sf = abs(fft2(Y_mesu(:,:,l))).^2; % image power spectrum
        fA = sum(Sf(:))/numel(Y_mesu(:,:,l)); % image average power
        R1 = abs(nA) / abs(fA);
        y1 = Y_mesu(:,:,l) +  noise;
        Y_rgb_noisy(:,:,l,k) = uint8(y1);

        if (algo==1)
            if(sigma_s(1) == 0)
                [wftempo(:,:,l)] = wienerFilter(PSF(:,:,l,k),double(Y_rgb_noisy(:,:,l,k)),deltaS,R1);
            else
                wftempo(:,:,l) = deconvwnr(double(Y_rgb_noisy(:,:,l,k)),PSF(:,:,l,k),R1);
            end
        elseif(algo==2)
            [wftempo(:,:,l)] = mat2gray(l2_TV(double(Y_rgb_noisy(:,:,l,k)),PSF(:,:,l,k)));
        end
    end

    wf(:,:,k,:)  = norm01(wftempo); %(wftempo ./max(wftempo(:)));
end
if(ploton==1)
    for k=1:K
        subplot(3,K+1,k+1),imagesc(uint8(Y_rgb_noisy(:,:,:,k))),title("Measurement z="+ num2str(z(k))+"[m]"),axis off,pbaspect([1 1 1])
    end
    subplot(3,K+1,2*(K+1)+1),imagesc(d),title('DOE'),pbaspect([1 1 1]),axis off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:K
    wftemp1 = uint8((squeeze(wf(:,:,k,:)))*255);
    p(k) = psnr(wftemp1,RGB_X);
    si(k) = ssim(wftemp1,RGB_X);
    p1 = p(k);
    s1 = si(k);
    if(ploton==1)
        subplot(3,K+1,K+2),plot(linspace(-Ns,Ns,N),d(Ns,:),color{3}),title("DOE profile"),axis([-Ns Ns 0 max(d(Ns,:))]);
        if(showMTF == 1)
            for l=1:L
                H(:,:,l)  = psf2otf(PSF(:,:,l,k));
                temp = abs(fftshift(H(:,:,l)));
                H1(:,:,l) = temp ./max(temp(:));
                subplot(3,K+1,(K+1)+k+1),plot(linspace(-1,1,Ns),H1(round(Ns/2),:,l),color{l}),title("MTF "+"z="+ num2str(z(k))+"[m]"),pbaspect([1 1 1]),axis on,hold on;
            end
            legend(wavelength(1)+'nm',wavelength(2)+'nm',wavelength(3)+'nm','FontSize',5,'Location','northwest')
        else
            for l=1:3
                h1(:,:,l)  = PSF(:,:,l,k);
            end
            h1 = norm01(h1);
            subplot(3,K+1,(K+1)+k+1),imagesc((h1).^gamma),title("PSF "+"z="+ num2str(z(k))+"[m]"),pbaspect([1 1 1]),axis off;
        end
        if(ploton==1)
            subplot(3,K+1,2*(K+1)+k+1),imagesc(wftemp1),title(algorithm(algo)+" z="+ num2str(z(k))+"[m]" +", PSNR="+ num2str(p(k))),pbaspect([1 1 1]),axis off;
            colormap('jet')
        end
    end
end
pm = mean(p);
sm = mean(si);
disp("PSNR= "+pm+" PSNR std= "+std(p));
disp("SSIM= "+sm+" SSIM std= "+std(si));

end
function x = norm01(x)
mini=min(x,[],"all");
maxi=max(x,[],"all");
x=(x-mini)/(maxi-mini);
end