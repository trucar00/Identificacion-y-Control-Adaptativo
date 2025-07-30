% Initialize
N = 100; % we run for 200 seconds

lambda = 0.99; % Forgetting factor (adjust between 0.98-1 for noise rejection)

% Initialize parameters
theta_hat = zeros(2, 1); % [a1, b1]
P = 1e5 * eye(2);  % Large initial covariance

alpha = -0.5; % alc1

ref = ones(N,1)*h_list(3);
u = zeros(N,1);
y = zeros(N,1);
y(1) = h_list(1);
u(1) = k_list(1);

theta_hist = zeros(2, N);

for k = 2:N
    
    time = (0:N-1)*0.1;
    u_timeseries = timeseries(u, time);
    % --- Step 1: Send Control Input to Simulink ---
    % Run Simulink for one time step
    simOut = sim('tp4_2_str.slx', 'StopTime', num2str(k));

    y_real = simOut.nivout.Data; % Output
    u_real = simOut.ctrout.Data; % Input
    y_sampled = y_real(1:10:end); % Sampled output 
    u_sampled = u_real(1:10:end); % Sampled input
    y_sampled = y_sampled - y_sampled(1);
    u_sampled = u_sampled - u_sampled(1);

    N_rls = length(y_sampled); % Number of samples

    for j = 2:N_rls
        x_k = [y_sampled(j-1); u_sampled(j-1)];  % Include past inputs
        
        y_hat = x_k' * theta_hat;
        
        K = (P * x_k) / (lambda + x_k' * P * x_k);  % Gain update with forgetting factor
        theta_hat = theta_hat + K * (y_sampled(j) - y_hat); % Update parameter estimate
        theta_hist(:, k) = theta_hat;
        P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k); % Covariance update
    end

    % Extract estimated parameters
    a1 = theta_hat(1);
    b1 = theta_hat(2);

    t0 = -(1+alpha)/b1;
    s0 = -(alpha-a1)/b1;

    y(k) = a1 * y(k-1) + b1 * u(k-1);

    u(k) = t0 * ref(k) - s0 * y(k);
end