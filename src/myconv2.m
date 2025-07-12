function C = myconv2(A, B, delta)
% function C = myconv2(A, B, delta)
M = size(A, 1);
N = size(A, 2);
%C = ift2(ft2(A, delta) .* fftshift(ft2(B, delta)), 1/(N*delta));
C = ift2((ft2(A, delta) .* ft2(B, delta)), 1/(N*M*delta));