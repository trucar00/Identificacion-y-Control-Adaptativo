% Initialize
N = 30; % we run for 200 seconds

n = 1; % Output order
m = 1; % Input order 
lambda = 0.99; % Forgetting factor (adjust between 0.98-1 for noise rejection)

% Initialize parameters
theta_hat = zeros(n + m, 1); % [a1, b1]
P = 10 * eye(n + m);  % Large initial covariance

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
    simOut = sim('tp4_2_str.slx', 'StopTime', num2str(k-1));
    
    % --- Step 2: Get Actual Output from Simulink ---
    y(k) = simOut.nivout.Data(end);

    % --- Step 3: Update Parameter Estimates (RLS) ---
    x_k = [y(k-1); u(k)]; % Use past values
    y_hat = x_k' * theta_hat; % Predicted output

    K = (P * x_k) / (lambda + x_k' * P * x_k); % RLS Gain
    theta_hat = theta_hat + K * (y(k) - y_hat); % Update parameters
    P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k); % Covariance update

    % Extract estimated parameters
    a1 = theta_hat(1);
    b1 = theta_hat(2);

    theta_hist(:, k) = theta_hat;

    % --- Step 4: Compute Updated Control Law ---
    t0 = -(1+alpha)/b1;
    s0 = -(alpha-a1)/b1;

    u(k) = t0 * ref(k) - s0 * y(k);
end