%
% Time series analysis
% Lund University
%
% Example code 18: Recursive least squares.
%
% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
clear; clc;
close all;

% Simulate an AR process. 
rng(1)                                          % Set the seed (just done for the lecture!)
extraN = 100;
N = 1000;
A = [ 1 -0.6 0.8 ];
%A = [ 1 -0.6 0.8 0.2 ];                        % Try changing the process.
e = randn( N+extraN, 1 );
y = filter( 1, A, e );     y = y(extraN:end);
p = length(A)-1;

% Plot realisation.
figure
plot(y)
title( sprintf('Realisation of an AR(%i)-process', p) )
ylabel('Amplitude')
xlabel('Time')
xlim([1 N])

figure
zplane(A)
title('Roots of the A(z) polynomial')


%% Estimate the parameters using LS (see example 5.6).
X = hankel( y(p:-1:1), y(1:N-p+1) ).';
theta = -inv( X.'*X )*X.'*y(p:N);
stdThLS = 2*std(e)*sqrt( diag(inv( X.'*X )) );
fprintf('Estimating the AR coefficients.\n')
fprintf('  True coefficients: %4.2f and %4.2f.\n', A(2), A(3) )
fprintf('  LS estimate:       %4.2f +/- %4.2f.\n', theta(2), 2*stdThLS(2) );
fprintf('  LS estimate:        %4.2f +/- %4.2f.\n', theta(1), 2*stdThLS(1) );


%% Form the RLS estimate (see section 8.3)
lambda = 1;                                     % This is the forgetting factor. Try, e.g., 0.95-0.999.
Pt1 = eye(p);                                   % Initial value of P_{t-1}.
th  = zeros(p,N);                               % Initial values of \theta_t.
for t=p+1:N
    xt   = -y(t-1:-1:t-p);                      % Form regressor (see also example 5.6).
    Pt1x = Pt1*xt;                              % Pre-compute to save calculations.
    nom  = lambda + xt.'*Pt1x;

    Kt  = Pt1x / nom;
    Pt1 = ( Pt1 - Pt1x*Pt1x.'/nom )/lambda;  

    th(:,t) = th(:,t-1) + Kt*( y(t) - xt.'*th(:,t-1) );
end

figure
plot( th.')
line([1 N], [A(2) A(2)], 'Color','red','LineStyle',':')
line([1 N], [A(3) A(3)], 'Color','blue','LineStyle',':')
legend('Estimate a_1', 'Estimate a_2', 'True a_1', 'True a_2','Location','best')
xlabel('Time')
ylabel('Estimated coefficients')
ylim([-1.5 1.5])
title( sprintf('RLS estimate with lambda = %3.2f.', lambda))
stdThRLS = sqrt( diag(Pt1) );
fprintf('Using a forgetting factor lambda = %4.3f.\n', lambda);
fprintf('  RLS estimate:      %4.2f +/- %4.2f.\n', th(1,N), 2*stdThRLS(1) );
fprintf('  RLS estimate        %4.2f +/- %4.2f.\n', th(2,N), 2*stdThRLS(2) );
