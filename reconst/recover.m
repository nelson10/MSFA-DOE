function [Xrec] = recover(Y,T)
K = size(Y,3);
[M,N,L] = size(T);
Xrec = zeros(M,N,L,K);
dynamicRange = 2^8-1;
for i=1:K
    [J] = unfold(Y(:,:,i),T);
    [Xrec(:,:,:,i)] = interpolation(J);
end
Xrec(Xrec<=0)=0;
Xrec(Xrec>=dynamicRange)=dynamicRange;
end