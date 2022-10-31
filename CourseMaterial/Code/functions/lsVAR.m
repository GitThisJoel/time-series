% FUNCTION [ thEst, sigEst, resEst ] = lsVAR( y, p )
%
% The function estimates the multivariate AR matrices for a VAR(p) process.
%

% Reference: 
%   "An Introduction to Time Series Modeling" by Andreas Jakobsson
%   Studentlitteratur, 2013
%
function [ thEst, sigEst, resEst ] = lsVAR( y, p )

% Arrange Y to be column vectors.
[m,N] = size(y);
if m>N
    y = y.';
    [m,N] = size(y);
end

% Build X and Y matrices.
Y = y( :,p+1:N )';
X = [];
for k=p+1:N   
    Xp = - y( :, k-1:-1:k-p);
    Xp = Xp(:)';
    X = [ X ; Xp ];
end
thEst  = inv( X'*X )*X'*Y;
resEst = Y-X*thEst;
sigEst = resEst' * resEst / (N-p);

thEst = thEst';