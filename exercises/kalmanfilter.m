function kalmanfilter(y)
    % Example of Kalman filter
    % Simulate N samples of a process to test your code .

    y = ? % Simulated data

    % Define the state space equations.

    A = [? ?; ? ?];
    Re = [? ?; ? ?]; % State covariance matrix
    Rw = ? % Observation variance

    % Set some initial values
    Rxx1 = ? * eye(2); % Initial state variance
    xtt1 = [? ?]'; % Initial statevalues

    % Vectors to store values in
    Xsave = zeros(2, N); % Store d states
    ehat = zeros(1, N); % Prediction residual
    yt1 = zeros(1, N); % One step prediction
    yt2 = zeros(1, N); % Two step prediction

    % The filter use data up to time t −1 to predict value at t,
    % then update using the prediction error. Why do we start
    % from t = 3? Why stop at N−2?
    for t = 3:N - 2
        Ct = [? ?]; % C { t | t −1}
        yhat(t) = ? % y { t | t −1}
        ehat(t) = ? % e t = y t - y { t | t −1}

        % Update
        Ryy = ? % Rˆ{ yy } { t | t −1}
        Kt = ? % Kt
        xtt = ? % x{ t | t }
        Rxx = ? % R{ xx } { t | t }

        % Predict the next state
        xtt1 = ? % x{ t+1| t }
        Rxx1 = ? % Rˆ{ xx }{ t+1 | t }

        % Form 2− step prediction. Ignore this part at first.
        Ct1 = [? ?]; % C{ t+1 | t }
        yt1(t +1) = ? % y{ t+1 | t } = C{ t+1 | t } x{ t | t }
        Ct2 = [? ?]; % C{ t+2 | t }
        yt2(t +2) = ? % y{ t+2 | t } = C{ t+2 | t } x{ t | t }

        % Store the state vector
        Xsave(:, t) = xtt;
    end

end
