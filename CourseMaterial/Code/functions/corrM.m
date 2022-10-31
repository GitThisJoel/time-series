% FUNCTION [Ry, rhoY] = corrM( Y, L )
%
% The function estimates the auto-covariance and auto-correlation
% estimates for lags 0 up to L for a multivariate process.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [Ry, rhoY] = corrM( Y, L )

[m,N] = size( Y );
if N<m
    Y = Y.';
    [m,N] = size( Y );
end
mY = mean(Y,2);

for k=0:L  
    Ry{k+1} = zeros(m,m);
    for p=k+1:N
        Ry{k+1} = Ry{k+1} + ( Y(:,p)-mY )*( Y(:,p-k)-mY )';
    end
    Ry{k+1} = Ry{k+1} / N;
end

D = diag( 1./(sqrt( diag(Ry{1}) )) );
for k=0:L  
    rhoY{k+1} = D*Ry{k+1}*D;
end
