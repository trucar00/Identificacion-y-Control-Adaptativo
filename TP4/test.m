% Initialize
N = 20; % we run for 200 seconds

lambda = 0.99; % Forgetting factor (adjust between 0.98-1 for noise rejection)

% Initialize parameters
theta_hat = zeros(2, 1); % [a1, b1]
P = 1e5 * eye(2);  % Large initial covariance

alpha = -0.5; % alc1

ref = ones(N,1)*h_list(3);

%theta_hist = zeros(2, N);

u = zeros(N,1);
u(1) = k_list(1);

%y_est = zeros(N,1);
%y_est(1) = h_list(1);


u_max = 1;  % Maximum control input
u_min = 0;  % Minimum control input

for k = 2:N
    
    time = (0:N-1);
    u_timeseries = timeseries(u, time);
    % --- Step 1: Send Control Input to Simulink ---
    % Run Simulink until k
    simOut = sim('tp4_2_str.slx', 'StopTime', num2str(k));

    y = simOut.nivout.Data(1:10:end);
    disp(y)
    u_rls = simOut.ctrout.Data(1:10:end); % Input
    y = y - y(1);
    u_rls = u_rls - u_rls(1);

    K = length(y);

    for j = 2:K
        x_k = [y(j-1); u_rls(j-1)];  % Include past inputs
        
        y_hat = x_k' * theta_hat;
        
        K = (P * x_k) / (lambda + x_k' * P * x_k);  % Gain update with forgetting factor
        theta_hat = theta_hat + K * (y(j) - y_hat); % Update parameter estimate
        %disp(theta_hat)
        P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k); % Covariance update
    end

    % Extract estimated parameters
    a1 = theta_hat(1);
    b1 = theta_hat(2);

    t0 = -(1+alpha)/b1;
    s0 = -(alpha-a1)/b1;

    %y_est(k) = a1 * y_est(k-1) + b1 * u(k-1);

    u(k) = t0 * ref(k) - s0 * y(k);
    u(k) = max(min(u(k), u_max), u_min);
    
end