function d = helicalDOE(N,doe_pitch,pars)
    height_max = 2.4e-6;
    f = pars(1);
    betha = pars(2);
    alpha = pars(3);
    L = pars(4); % number of spirals
    
    V=0.5;%%Visibility controller
    lambda=0.550*1e-6;%Define wavelength
    del=doe_pitch; % doe pitch
    x=-N/2:N/2-1;
    y=-N/2:N/2-1;
    [X,Y]=meshgrid(x*del,y*del);
    R=sqrt(X.^2+Y.^2);

    d = V*exp(alpha*1i*L*(atan2(X,Y))  + betha*1i*((2*pi)/lambda)*(f-sqrt(f*f-R.*R)));
    d = mat2gray(imag(d));
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    N1 = 128;
    d = mod(d*N1 , N1);
    d = floor(d);
    d = mat2gray(d) .* height_max;
    d = d .*A_st;
end
