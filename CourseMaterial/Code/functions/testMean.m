% FUNCTION [ rejectMean, tRatio, tLimit] = testMean( data, mean0, signLvl )
%
% The function tests if the mean of the data can be deemed to be that of a 
% process with mean-value mean0 (default zero), with significance signLvl 
% (default 0.05), returning the decision, as well as the test ratio and the  
% decision limit (for univariate data this is t-distributed, for multi-variate, 
% it is Hotelling-T2 distributed).
%
% The function assumes that the data is  Gaussian distributed, and works 
% for both univariate and multivariate data sets, with the latter assumed 
% to consist of N row vectors, formed as an N x m matrix. 
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ rejectMean, tRatio, tLimit] = testMean( data, mean0, signLvl )


[N,m] = size( data );
if (N==1) && (m>1)     % Univariate data given as row vector.
    N = m;
    m = 1;
end
if nargin<2    
    mean0 = zeros( 1, m );
end
if nargin<3
    signLvl = 0.05;
end

if m==1                % Form squared t-ratio and the rejection limit.
    tLimit = tinv( 1-signLvl/2, N-1 ).^2;
    tRatio = N*( mean(data)-mean0 )^2 / var(data);
    
else                    % Form Hotelling T2 and the rejection limit.
    R0 = covMvect( data, mean0 );
    R1 = covMvect( data, mean(data) );
    tLimit = (N-1)*m/(N-m) * finv( 1-signLvl, m, N-m );
    tRatio = (N-1)*det( R0 ) / det( R1 ) - (N-1);
end
rejectMean = tRatio > tLimit;
