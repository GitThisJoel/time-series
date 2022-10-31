% FUNCTION [ isNormal, P ] = dagostinoK2test( data, signLvl )
%
% The function computes the D'Agostino-Pearson's K2 test of normality with
% significance signLvl (default 0.05). 
%

% The function is a slightly modified version of the function provided by 
% A. Trujillo-Ortiz and R. Hernandez-Walls, which is available as:
%
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). DagosPtest: D'Agostino-Pearson's K2 test for 
%    assessing normality of data using skewness and kurtosis. A MATLAB file. [WWW document]. URL 
%    http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=3954&objectType=FILE
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [ isNormal, P ] = dagostinoK2test( data, signLvl )

if nargin<2
    signLvl = 0.05;
end
if (signLvl <= 0 || signLvl >= 1)
   error('dagostinoK2test: significance level must be between 0 and 1\n');
end

n     = length(data);
data  = sort( data(:) );

[c,v] = hist( data, data );
nc    = find( c~=0 );
c     = [ v(nc) c(nc)' ];
x     = c(:,1);
f     = c(:,2);
s1    = f'*x;
s2    = f'*x.^2;
s3    = f'*x.^3;
s4    = f'*x.^4;
SS    = s2-(s1^2/n);
v     = SS/(n-1);
k3    = ((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2));
g1    = k3/sqrt(v^3);
k4    = ((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(SS^2))/((n-2)*(n-3)));
g2    = k4/v^2;
eg1   = ((n-2)*g1)/sqrt(n*(n-1));                         % This is a measure of the skewness
eg2   = ((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1)); % This is a measure of the kurtosis
A     = eg1*sqrt(((n+1)*(n+3))/(6*(n-2)));
B     = (3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9));
C     = sqrt(2*(B-1))-1;
D     = sqrt(C);
E     = 1/sqrt(log(D));
F     = A/sqrt(2/(C-1));
Zg1   = E*log(F+sqrt(F^2+1));
G     = (24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5));
H     = ((n-2)*(n-3)*g2)/((n+1)*(n-1)*sqrt(G));
J     = ((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))));
K     = 6+((8/J)*((2/J)+sqrt(1+(4/J^2))));
L     = (1-(2/K))/(1+H*sqrt(2/(K-4)));
Zg2   = (1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K));

K2    = Zg1^2 + Zg2^2;
P     = 1-chi2cdf(K2,2);

isNormal = ( P >= signLvl );
