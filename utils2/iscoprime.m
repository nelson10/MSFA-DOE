function [iscp, ispcp, pidx, pgcd] = iscoprime(x)
% ISCOPRIME     Check coprime relation
%   ISCP = iscoprime(X) returns the logical value ISCP as true if the
%   elements in X are coprime, and as false if the elements in X are not
%   coprime. The elements in X are coprime when the greatest common divisor
%   of all elements in X is 1. X is a M-element row vector of integers.
%
%   [ISCP,ISPCP,PIDX,PGCD] = iscoprime(X) also returns ISPCP to indicate if
%   the elements in X are pairwise coprime. When ISPCP is true, any two
%   elements in X are coprime. PIDX is a 2-row matrix whose number of
%   columns equals to M choose 2. Each column of PIDX specifies the
%   indices of a pair of elements in X. PGCD is a row vector with the same
%   number of columns as PIDX. Each element in PGCD is the greatest common
%   divisor of the two elements in X identified by the indices in the
%   corresponding column of PIDX.
%
%   % Example:
%   %   Check if the vector [21 36 49] is coprime. Is it pairwise coprime?
%
%   [iscp,ispcp] = iscoprime([21 36 49])
%
%   See also phased, crt.

%   Copyright 2020 The MathWorks, Inc.

%#codegen

% validation
validateattributes(x,{'double','single'},{'nonnegative','row',...
    'nonnan','finite','nonempty','integer'},'iscoprime','X');
cond = numel(x)<2;
coder.internal.errorIf(cond,'phased:phased:tooFewColumns','X','2','IfNotConst','CheckAtRunTime');

pidx = nchoosek(1:numel(x),2).';
pgcd = gcd(x(pidx(1,:)),x(pidx(2,:)));
iscp = any(pgcd==1);
ispcp = all(pgcd==1);


% LocalWords:  coprime ISCP ISPCP PGCD iscp ispcp crt PIDX
