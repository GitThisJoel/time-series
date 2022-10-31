% FUNCTION [ An,Bn,Pf,Pb ] = wwra( Rfb )
%
% The function implements the Whittle-Wiggins-Robinson Algorithm
% (WWRA). The function takes M1+1 Block matrices as input (which
% represents the full Block-Toeplitz matrix) and gives the forward
% and backward prediction matrices as output.
%
% Reference : 
%   Söderström & Stoica, "System identification", pp. 310-319. 
%
% Input :  Rfb - M1+1 (M2+1)x(M2+1) block matrices. 
%
% Output : An  - M1+1 (M2+1)x(M2+1) forward prediction matrices.
%          Bn  - M1+1 (M2+1)x(M2+1) backward prediction matrices.
%          Pf  - (M2+1)x(M2+1) forward prediction covariance matrix.
%          Pb  - (M2+1)x(M2+1) backward prediction covariance matrix.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ An,Bn,Pf,Pb ] = wwra( Rfb, M1, M2 )

% Initialization.
  An = zeros(M2+1,M2+1,M1+1);
  Bn = An;  
  Ao = An; 
  Bo = An;
  I  = eye(M2+1);
  
% WWRA initialization.
  An(:,:,1) = I;
  Bn(:,:,1) = I;
  Pf = Rfb(:,:,1);
  Pb = Pf;

% WWRA order update.
  for n=0:M1-1
    Ao = An;
    Bo = Bn;
    
    Pn = Rfb(:,:,n+2);
    for k=1:n
      Pn = Pn + Ao(:,:,k+1)*Rfb(:,:,n-k+2);
    end

    An(:,:,n+2) = - Pn*inv( Pb );
    Bn(:,:,n+2) = - Pn'*inv( Pf );
    for k=1:n
      An(:,:,k+1) = Ao(:,:,k+1) + An(:,:,n+2)*Bo(:,:,n-k+2);
      Bn(:,:,k+1) = Bo(:,:,k+1) + Bn(:,:,n+2)*Ao(:,:,n-k+2);
    end
    
    Pf = Pf + An(:,:,n+2)*Pn';
    Pb = Pb + Bn(:,:,n+2)*Pn;
  end 


