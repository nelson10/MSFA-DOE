function showPlot(X,wf,Y_rgb,d,PSF,j,showMTF,ploton,algo,img)
algorithm = ["Lucy","l2-TV"];
a = 1e0; % 1e0 [m]  1e3 [mm]
z = [0.5 1.0 2.0 inf].*a; %%
K =size(wf,4);
N =size(d,2);
Ns = round(N/2);
L =size(X,3);
colors = colormap(jet(L));
dynamicRange = 2^16-1;

wln = linspace(4.00e-7,7.00e-7,31); %linspace(4.82e-7,6.11e-7,31);
temp = wln*1e9;
wavelength = string(temp);

wf = mat2gray(double(wf)).*max(img(:));
X = mat2gray(double(X)).*max(img(:));
for k=1:4
    [p(1,k),s(1,k),sam(1,k)] = metrics2(X,wf(:,:,:,k));
end

if(ploton==1)
    rgb = RGB_test(X);
    subplot(3,K+1,1),imagesc(rgb(1:Ns,1:Ns,[3 2 1])),title("Groundtruth Image "+num2str(j)),pbaspect([1 1 1]),axis off;
end
if(ploton==1)
    for k=1:K
        rgb = RGB_test(Y_rgb(:,:,:,k));
        subplot(3,K+1,k+1),imagesc(rgb(:,:,[3 2 1])),title("Measurement z="+ num2str(z(k))+"[m]"),axis off,pbaspect([1 1 1])
    end
    subplot(3,K+1,2*(K+1)+1),imagesc(d),title('DOE'),pbaspect([1 1 1]),axis off;
end
for k=1:K
    wftemp1 = wf(:,:,:,k);
    % p(k) = psnr(wftemp1,X);
    % si(k) = ssim(wftemp1,X);
    % p1 = p(k);
    % s1 = si(k);
    if(ploton==1)
        subplot(3,K+1,K+2),plot(linspace(-Ns,Ns,N),d(Ns,:)),title("DOE profile"),axis([-Ns Ns 0 max(d(Ns,:))]);
        if(showMTF == 1)
            for l=1:L
                H(:,:,l)  = psf2otf(PSF(:,:,l,k));
                temp = abs(fftshift(H(:,:,l)));
                H1(:,:,l) = temp ./max(temp(:));
                subplot(3,K+1,(K+1)+k+1),plot(linspace(-1,1,Ns),H1(round(Ns/2),:,l),'Color',colors(l,:)),title("MTF "+"z="+ num2str(z(k))+"[m]"),pbaspect([1 1 1]),axis on,hold on;
                %legend(wavelength(l)+'nm','FontSize',5,'Location','northwest')
            end
        else
            for l=1:L
                h1(:,:,l)  = PSF(:,:,l,k);
            end
            h1 = mat2gray(h1).*dynamicRange;
            rgb = RGB_test(h1);
            subplot(3,K+1,(K+1)+k+1),imagesc((rgb(:,:,[3 2 1]).^0.25)),title("PSF "+"z="+ num2str(z(k))+"[m]"),pbaspect([1 1 1]),axis off;
            %pause(0.1);
            colormap('jet')

        end
        if(ploton==1)
            rgb = RGB_test(wftemp1);
            subplot(3,K+1,2*(K+1)+k+1),imagesc(rgb(:,:,[3 2 1])),title(algorithm(algo)+" z="+ num2str(z(k))+"[m]" +", PSNR="+ num2str(p(k))),pbaspect([1 1 1]),axis off;
            colormap('jet')
        end
    end
end
NL = 6;
idx = round(linspace(1,L,NL));
PSF = mat2gray(PSF).*dynamicRange;
figure(2)
colormap('jet')
for l=1:NL
    for k=1:K
        subplot(NL,K,(K)*(l-1)+k),imagesc(PSF(:,:,idx(l),k).^0.25),title('\lambda= '+wavelength(idx(l)) +' z='+num2str(z(k))),pbaspect([1 1 1]),axis off;
    end
end
end