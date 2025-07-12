% Created by: Ph.D Nelson Eduardo Diaz Diaz
% Post-doctorado Pontícia Universidad Católica de Valparaíso (PUCV)
% Date 2 February 2022

% Comparison of minimum distance
% Solution using the Discrete Sphere Packing based on 3D N^2 Queens
% Approach

% Find Optimal values of a and b

function [a,b,diameter,G1]=DDDRSNNPLattice(N,NF)
if(NF >= 16)
    M = round(NF/2);
    id = 1:M;
    t = zeros(M,1);
    for m = 1:M
        t(m) = iscoprime([m N]);
    end
    idx =id(t==1);
else
    M = NF;
    idx = 1:NF;
end
K = length(idx);

distance = zeros(M,M);
x = ones(1,NF)';
y = (1:NF)';
I = kron(x',y);
J = kron(x,y');

densityM = zeros(M,M);
distance1 = zeros(M,M);
distance2 = zeros(M,M);
n = 3;

for i=1:K
    %disp(num2str(i)+" out of "+num2str(K));
    for j=i:K

        si = idx(i);
        sj = idx(j);
        G = mod(I.*si + J.*sj,NF)+1;

        %
        [flag] = VerifyRollingShutter(G,NF);
        flag=1;
        if(flag == 1)
            %imagesc(G==1)
            %pause(2)
            deno = (NF.^2);
            nume =(NF+1).^3;
            [r,c,z] = find(G(1:NF,1:NF));
            X = [r c z]'; % generator matrix
            %if(NF<16)
            [r1,c1,z1] = find(G(1:NF,1:NF)==1);
            X1 = [r1 c1 z1]'; % generator matrix
            D = pdist(X1');
            NAlpha1 = min(D);
            distance1(si,sj) = NAlpha1;
            %end
            %M1 = orth(X); % Compute an orthonormal basis
            [U,S,V] = svd(X,'econ');
            A = U'*U;
            volC = sqrt(det(A)); % Volume of a fundamental cell
            %D = pdist(X');
            if(NF < 128)
                Q = round(U*S*V(1:round(length(X)),:)');
            else
                Q = round(U*S*V(1:round(0.050*length(X)),:)');
            end
            D = pdist(Q');
            NAlpha = min(D); % minimal value between <v,v> for all v \in \Alpha
            r1 = NAlpha/2;
            num = ((pi).^(n/2)) .* (r1.^n);
            den = gamma(n/2 + 1);
            volS= num / den; % volume of an n-ball, https://en.wikipedia.org/wiki/Volume_of_an_n-ball
            density = (volS/ volC)*(deno/nume);
            %density = ComputeDensity(NAlpha,NF);
            densityM(si,sj) = density;
            distance2(si,sj) = NAlpha;
            %[distance(si,sj)]=distan3(G,NF);

        end
    end
end

ma = max(densityM(:));

d = 1*distance1 + distance2;
maxi = max(d(:));
diameter = max(distance2(:));
%disp("diameter= "+num2str(diameter));
[b,a,~]= find(maxi==d);
G1 = mod(I.*a(1) + J.*b(1),NF)+1;
A = G1;
p = ceil(N/NF);
B = ones(p,p);
G1 = kron(B,A);
G1 = G1(1:N,1:N);
end
function [R] = bravaisLattice(G,N)
R = [];
for i=0:N
    for j=0:N
        for k=0:N
            z = [i,j,k];
            x = z*G;
            R = [R;x];
        end
    end
end

end