% FUNCTION [ S, Ry ] = pacfM( Y, L )
%
% The function estimates the partial auto-correlation
% estimates for lags 0 up to L for a multivariate process.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ S, Ry ] = pacfM( Y, L )

% Compute correlations.
Ry = corrM( Y, L );

% Find the partial lag correlation matrices. 
% First iteration, s=1
Vu{1}  = Ry{1};
Vv{1}  = Vu{1};
Vvu{1} = Ry{2};
tmp    = inv( Ry{1} );
alp{1,1} = Ry{2}' * tmp;
bet{1,1} = Ry{2}  * tmp;

Dv = diag( 1./sqrt( diag( Vv{1} )) );
Du = diag( 1./sqrt( diag( Vu{1} )) );
S{1} = Dv * Vvu{1}' * Du;

% For s>1,
for s=2:L-1
    Vu{s} = Ry{1};
    Vv{s} = Ry{1};
    Vvu{s} = Ry{s+1};
    for k=1:s-1
        Vu{s} = Vu{s} - alp{s-1,k}*Ry{k+1};
        Vv{s} = Vv{s} - bet{s-1,k}*Ry{k+1}';
        Vvu{s} = Vvu{s} - Ry{s-k+1}*alp{s-1,k}';      
    end

    alp{s,s} = Vvu{s}' * inv( Vv{s} );
    bet{s,s} = Vvu{s}*inv( Vu{s} );
 
    for k=1:s-1
        alp{s,k} = alp{s-1,k} - alp{s,s}*bet{s-1,s-k};
        bet{s,k} = bet{s-1,k} - bet{s,s}*alp{s-1,s-k};
    end 
        
    Dv = diag( 1./sqrt( diag( Vv{s} )) );
    Du = diag( 1./sqrt( diag( Vu{s} )) );
    
    S{s} = Dv * Vvu{s} * Du;
end
