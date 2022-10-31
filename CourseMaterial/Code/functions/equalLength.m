% FUNCTION varargout = equalLength(varargin)
%
% The function extends all entered polynomials such that the returned
% versions all have the same length.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function varargout = equalLength(varargin)

varLength = zeros(nargin,1);
for i=1:nargin
  varLength(i)=length(varargin{i});
end
maxLength = max(varLength);
for i=1:nargin
    out = zeros(1,maxLength);
    out(1:varLength(i)) = varargin{i};
    varargin{i} = out;
    
end
varargout =varargin;
