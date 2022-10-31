%
% Time series analysis
% Lund University
%
% Example code 19: State space model of an ARMA process.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate an ARMA process. Leave the intial values to see how the
% different representation affect these samples.
N  = 200;
C0 = [ 1 0.6 0.5 -0.8 ];
A0 = [ 1 -0.6 0.8 ];
e = randn( N, 1 );
y = filter( C0, A0, e );                        % This is the transfer function representation.

plot( y )
title('Transfer function representation')
xlabel('Time')


%% Form the state space representations of the ARMA.
% Note how these representations are equivalent (except in the intial "ring-down").
p = length(A0)-1;
q = length(C0)-1;
d = max(p, q+1);
fprintf('Number of states required is max(p, q+1) = %i.\n', d)

% Controllable canonical form.
A  = [ -A0(2) -A0(3) 0 0 ; 1 0 0 0; 0 1 0 0 ; 0 0 1 0];
B  = [ 1 ; 0 ; 0 ; 0];
C  = C0;
xt = [0 ; 0; 0; 0];
yc = zeros(N,1);
for t=1:N
    xt    = A*xt + B*e(t);
    yc(t) = C*xt;
end
fprintf('The A matrix for a controllable canonical form is:\n')
A


%% Observable canonical form.
A  = [ -A0(2) 1 0 0 ; -A0(3) 0 1 0; 0 0 0 1 ; 0 0 0 0 ];
B  = C0';
C  = [1 0 0 0];
xt = [0 ; 0; 0; 0];
yo = zeros(N,1);
for t=1:N
    xt    = A*xt + B*e(t);
    yo(t) = C*xt;
end
fprintf('The A matrix for an observable canonical form is:\n')
A


%% Random transformation of the state space form.
T  = randn(d,d);    % Create a random transformation matrix. 
iT = inv(T);
A  = T*A*iT;        % Transform the state matrices.
B  = T*B;
C  = C*iT;
yr = zeros(N,1);
for t=1:N
    xt    = A*xt + B*e(t);
    yr(t) = C*xt;
end
fprintf('The A matrix using a random transformation is:\n')
A


%% Compare realisations. 
figure
plot( [y yc yo yr] )
title('State space vs transfer representation')
legend('y_t', 'Controllable canonical form', 'Observable canonical form', 'Random form','Location','best')
xlabel('Time')
