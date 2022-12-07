function [xtt, ehat, Xsave] = kalmanfilter(y, A, C, N, Rw, Re, xtt1, Rxx1, ut)
    % Fixa sen: Ct
    
    % Define the state space equations.

    % Vectors to store values in
    Xsave = zeros(2, N); % Stored states
    yhat = zeros(1, N); %  Prediction residual
    ehat = zeros(1, N); %  Prediction residual
    yt1 = zeros(1, N); %   One step prediction
    yt2 = zeros(1, N); %   Two step prediction
    
    % The filter use data up to time t−1 to predict value at t,
    % then update using the prediction error. Why do we start
    % from t = 3? Why stop at N−2?
    for t = 3:N - 2
        Ct = [1 ut(t)]; % C{ t | t-1}
        yhat(t) = Ct * xtt1; % y{ t | t-1}
        ehat(t) = y(t) - yhat(t); % et = yt - y{ t | t-1 }
    
        % Update
        Ryy = Ct * Rxx1 * Ct' + Rw; % Rˆ{ yy }{ t | t-1 }
        Kt = Rxx1 * Ct' / Ryy; %  Kt
        xtt = xtt1 + Kt * (ehat(t)); % x{ t | t }
        Rxx = Rxx1 - Kt * Ryy * Kt'; % R{ xx }{ t | t }
    
        % Predict the next state
        xtt1 = A * xtt; % x{ t+1 | t }
        Rxx1 = A * Rxx * A' + Re; % Rˆ{ xx }{ t+1 | t }
    
        % Form 2− step prediction. Ignore this part at first.
        Ct1 = [-y(t) -y(t - 1)]; %   C{ t+1 | t }
        yt1(t + 1) = Ct1 * xtt; % y{ t+1 | t } = C{ t+1 | t } x{ t | t }
        Ct2 = [-yt1(t + 1) -y(t)]; %   C{ t+2 | t }
        yt2(t + 2) = Ct2 * xtt; % y{ t+2 | t } = C{ t+2 | t } x{ t | t }
    
        % Store the state vector
        Xsave(:, t) = xtt;
    end

end
