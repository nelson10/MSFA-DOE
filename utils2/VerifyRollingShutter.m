function [flag] = VerifyRollingShutter(G,N)
%Vflag = zeros(L,1);
temp0= zeros(N,1);
temp1= zeros(N,1);
j=1;
flag = 0;
s = sum(1:N);
while(j <= N)
    temp0(j) = sum(G(:,j))==s;
    temp1(j) = sum(G(j,:))==s;
    if(temp0(j)*temp1(j)==1)
        flag = 1;
    else
        flag = 0;
        break;
    end
    j=j+1;
end
end