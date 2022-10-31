% FUNCTION [ R, meanD ] = covMvect( data, [meanD], [fb] )
%
% The function computes the unbiased covariance matrix estimate of the data, 
% which is assumed to consist of N row vectors, formed as an N x m matrix. 
% The data is assumed to have a 1 x m mean vector meanD; if not given, it is 
% estimated. If fb is set, forward-backward averaging is used.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ R, meanD ] = covMvect( data, meanD, fb )

[N,m] = size( data );
if nargin<2                    % If not specified, estimate the mean.
    meanD = mean(data);
else
    [n1,m1] = size( meanD );    % Check dimensionality.
    if (n1==m) && (m1==1)       % If flipped, transpose to row vector.
        meanD = meanD.';
        m1 = n1;
        n1 = 1;
    end
    if m1~=m                   % Dimensionality error.
        error('covMvect: incompatible dimensions for mean vector and data.');
    end
end
if nargin<3
    fb = 0;
end

data =  data - ones(N,1)*meanD; % Compensate for mean value.
R    = data'*data / (N-1);      % Form forward-only covariance matrix.
if fb 
    H = fliplr(eye(m));         % Exchange matrix
    R = (R + H*R.'*H)/2;        % Form forward-backward estimate.
end
