function [Xrec] = interpolation(J)
[M,N,L]= size(J);
Xrec = zeros(M,N,L);
md = 5;
for j=1:L
    X = J(:,:,j);
    [x,y,v] = find(X);
    [xq,yq] = find(X==0);
    %x1 = interp2(x,y,v,xq,yq);
    F = scatteredInterpolant(x,y,v,'natural');
    %F1 = interp2(x,y,v);
    %X(X==0) = F1(xq,yq);
    X(X==0) = F(xq,yq);
    Xrec(:,:,j) = X; 
    %Xrec(:,:,j) = medfilt2(Xrec(:,:,j),[md md]);
    Xrec(:,:,j) = imsharpen(Xrec(:,:,j));
end
end