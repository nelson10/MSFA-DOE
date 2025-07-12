function [p,s,sam] = metrics2(data,Xrec)
[M,N,L] = size(Xrec);

p = zeros(L,1);
s = zeros(L,1);
rm = zeros(L,1);
for i=1:L
    ref = data((1+6):(N-6),(1+6):(N-6),i);
    A = Xrec((1+6):(N-6),(1+6):(N-6),i);
    p(i) = psnr(A,ref);
    s(i) = ssim(A,ref);
    v1 = double(A(:));
    v2 = double(ref(:));
    rm(i) = sqrt(immse(v1(:),v2(:)));
end

sa = zeros(M,N);
for n=1:N
    for m=1:M
        v1 = Xrec(m,n,:);
        v2 = data(m,n,:);
        v1 = double(v1(:));
        v2 = double(v2(:));
        sa(m,n) = real(SpectralAngleMapper(v1,v2));
        %tp(m,n) = isreal(sa(m,n));
    end
end

p = mean(p); % psnr
s= mean(s); % ssim
%r = mean(rm); % rmse
sam = mean(sa(:)); % sam Measure spectral similarity using spectral angle mapper
%+" RMSE "+num2str(r)
disp("PSNR "+ num2str(p)+ " SSIM "+num2str(s)+" SAM "+num2str(sam));
disp("---------------------------------------------------------------------------------------------------------")
end