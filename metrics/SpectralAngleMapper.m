%% Measure spectral similarity using spectral angle mapper

function [ang] = SpectralAngleMapper(v1,v2)
    v1 = v1 + eps;
    v2 = v2 + eps;    
    den = sqrt(sum(v1.^2)) * sqrt(sum(v2.^2));
    num = dot(v1,v2);
    ang = acos((num)/(den));
%     if(isreal(num)==0)
%     v1
%     v2
%     end
    %ang = ang/pi;
end