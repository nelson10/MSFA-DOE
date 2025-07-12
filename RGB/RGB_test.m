function [RGB]= RGB_test(dataset)

load illum_4000.mat; % evening sunlight with CCT 4000 K 
load illum_6500.mat; % daylight with CCT 6500 K
load('illum_25000.mat') % the spectra of blue skylight with correlated colour temperature (CCT) 25000 K
load xyzbar.mat; % CIE 1931 colour-matching functions

cube = dataset;
L2 = size(cube,3);

idx = round(linspace(1,L2,L2));
Io = mat2gray(cube(:,:,idx));
MS = Io;
%load(database); % Data
%hyperimg = reflectances(1:512,1:512,:); %Ori %bestRes
hyperimg = double(MS);
[N1,N2,L]=size(hyperimg);
reflectances=hyperimg;
l=round(linspace(1,33,L));

radiances_6500 = zeros(size(reflectances));  % initialize array
for j=1:L
    temp = reflectances(:,:,j);
    reflectances(:,:,j) = reflectances(:,:,j)/max(temp(:));
    radiances_6500(:,:,j) = reflectances(:,:,j);%*illum_6500(l(j)); % daylight with CCT 6500 K
end
[r, c, w] = size(radiances_6500);
radiances_6500 = reshape(radiances_6500, r*c, w);

XYZ = (xyzbar(l(:),:)'*radiances_6500')';
XYZ = reshape(XYZ, r, c, 3);
XYZ = max(XYZ, 0);
XYZ = XYZ/max(XYZ(:));

RGB = XYZ2sRGB_exgamma(XYZ);
RGB = max(RGB, 0);
RGB =RGB ./max(RGB(:));
RGB = min(RGB, 1);
gamma = 0.5;
pos   = [196 225;75 282;231 332];%;7 102;79 116]; % position bear
%pos   = [105 53]%; 9 15; 22 120];
color = {'red'};
%RGB = insertMarker(RGB,pos,'o','color',color,'size',5);
%figure(2); imshow((RGB).^gamma, 'Border','tight');
%imwrite(RGB.^gamma, "RGB"+num2str(k)+".png");
RGB = RGB.^gamma;
end