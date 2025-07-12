function [J] = unfold(Y,T)
[M,N,L] = size(T);
J = zeros(M,N,L);
for l=1:L
    J(:,:,l)= Y.*T(:,:,l);
end
end