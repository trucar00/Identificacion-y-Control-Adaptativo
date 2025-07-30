N = 200; % Run for 200 seconds
lambda = 0.99; % Forgetting factor (adjust for noise rejection)

% Initialize system parameter estimates
n = 1; % Order of output
m = 1; % Order of input
theta_hat = zeros(n + m, 1); % [a1, b1]
P = 1e5 * eye(n + m);  % Large initial covariance

% Initialize STR parameters
ref = ones(N,1) * h_list(3);  % Reference level
u = zeros(N,1);
y = zeros(N,1);
y(1) = h_list(1); % Initial height
u(1) = k_list(1); % Initial valve opening

u_max = 1;  
u_min = 0;  

% ======= MAIN LOOP =======
for k = 2:N
    
    % === 1. Create Timeseries for Simulink Input ===
    time = (0:k-1) * 10;  % Adjust the time vector dynamically
    u_timeseries = timeseries(u(1:k), time);  % Send only up to current step

    % === 2. Run Simulink for Current Time Step ===
    simOut = sim('tp4_2_str.slx', 'StopTime', num2str(k * 10));  % Adjust time scale

    % Retrieve measured data
    y_real = simOut.nivout.Data;
    u_rls = simOut.ctrout.Data;

    % Downsample and normalize
    y_measured = y_real(1:10:end);
    u_rls = u_rls(1:10:end);
    y_measured = y_measured - y_measured(1);
    
    % Use the latest measurement
    y(k) = y_measured(end);

    % === 3. Recursive Least Squares (RLS) Estimation ===
    K = length(y_measured);
    for j = 2:K
        x_k = [y_measured(j-1); u_rls(j-1)];
        
        y_hat = x_k' * theta_hat;
        
        K_rls = (P * x_k) / (lambda + x_k' * P * x_k);
        theta_hat = theta_hat + K_rls * (y_measured(j) - y_hat);
        P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k);
    end

    % Extract estimated parameters
    a1 = theta_hat(1);
    b1 = theta_hat(2);

    % Regularize b1 to avoid division by zero
    eps = 1e-3;
    b1 = b1 + eps;

    % === 4. Self-Tuning Regulator (STR) Control Law ===
    alpha = -0.5;
    t0 = -(1+alpha) / b1;
    s0 = -(alpha-a1) / b1;

    u(k) = t0 * ref(k) - s0 * y(k);
    disp(u)
    % Apply saturation
    u(k) = max(min(u(k), u_max), u_min);

    % Debugging: Print estimated parameters and control input
    %fprintf('Step %d: a1 = %.4f, b1 = %.4f, u = %.4f\n', k, a1, b1, u(k));
    
end