function [conditionNumber] = ComputeConditionNumber(T)
[N,M,L] = size(T);
kappa = zeros(N*M,L);
for j=1:L
    temp = T(:,:,j);
    kappa(:,j) = temp(:);
end

conditionNumber = cond(kappa'*kappa);
%disp("Codition Number "+ conditionNumber)
end