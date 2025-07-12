function [density,diameter]= ComputeDensityIrregularSP(T)
NF = size(T,3);
NF2 = size(T,2);
X = [];
for i=1:NF
    tp = T(1:NF2,1:NF2,i)*i;
    [x,y,z] = find(tp);
    temp =  [x,y,z];
    X = [X; temp];
end
%[r,c,z1] = find(G);
%X = [r c z1]'; % generator matrix
X = X';
[U,S,V] = svd(X,'econ');
if(NF <= 48)
    Q = round(U*S*V(1:round(length(X)),:)');
else
    Q = round(U*S*V(1:round(0.050*length(X)),:)');
end
D = pdist(Q');
D1 = squareform(D);
R = zeros(length(X),1);
for i=1:length(X)
     temp = D1(i,:);
     [x2,y2,z2]=find(temp);
     R(i)= min(z2)/2;
end

cube = (NF2+1).^3;
%dist = dist./2;
sphere = 0;
for i =1:length(X)
    sphere = sphere + (4/3).*pi.*(R(i).^3);
end
density = (sphere./cube);
diameter = mean(R(:))*2;
disp("Sphere Packing Density "+ density + " diameter "+ diameter);

% nb =NF.^2;
% dis = 2*(((NF+1).^3)./(4*nb*sqrt(2))).^(1/3);
% dis = dis./2;
% sphere = (NF.^2).*(4/3).*pi.*(dis.^3);
% density1 = (sphere./cube);
% disp("Optimal Sphere Packing Density "+ density1)
end