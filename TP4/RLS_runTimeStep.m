N = 200; % we run for 200 seconds

lambda = 0.99; % Forgetting factor (adjust between 0.98-1 for noise rejection)

% Initialize parameters
theta_hat = zeros(2, 1); % [a1, b1]
P = 1e5 * eye(2);  % Large initial covariance

alpha = -0.5; % alc1

ref = ones(N,1)*h_list(3);

theta_hist = zeros(2, N);

u = ones(N,1)*k_list(3);
u(1) = k_list(1);

for k = 2:N
    
    time = (0:N-1);
    u_timeseries = timeseries(u, time);
    % --- Step 1: Send Control Input to Simulink ---
    % Run Simulink until k
    simOut = sim('tp4_2_str.slx', 'StopTime', num2str(k));
    
    y = simOut.nivout.Data(1:10:end);
    u_rls = simOut.ctrout.Data(1:10:end); % Input
    
    %y = y - y(1);
    %u_rls = u_rls - u_rls(1);

    NN = length(y);

    for j = 2:NN
        x_k = [y(j-1); u_rls(j-1)];  % Include past inputs
        
        y_hat = x_k' * theta_hat;
        
        K = (P * x_k) / (lambda + x_k' * P * x_k);  % Gain update with forgetting factor
        theta_hat = theta_hat + K * (y(j) - y_hat); % Update parameter estimate
        
        P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k); % Covariance update
    end
    theta_hist(:, k) = theta_hat;
    % Extract estimated parameters
    a1 = theta_hat(1);
    b1 = theta_hat(2);

end

% Predicted output using estimated parameters
yh = zeros(NN, 1);
for i = 1:NN-1
    yh(i+1) = a1 * yh(i) + b1 * u_rls(i);
end

y_real = out.nivout.Data;
t_real = out.tout;
t = t_real(1:10:end);
figure;
plot(t_real, y_real, 'r', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(t, yh, 'b--', 'DisplayName', 'Identificación', LineWidth=1.8);
legend;
title('Comparación de salida, método recursivo');
xlabel('Tiempo');
ylabel('Salida');
fontsize(16, "points");