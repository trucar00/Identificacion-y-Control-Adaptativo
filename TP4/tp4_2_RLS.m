% Recursive Least Squares (RLS) with Forgetting Factor

y = out.nivout.Data(1:10:end); % Output (with noise)
u = out.ctrout.Data(1:10:end); % Input

% Subtract initial values
y = y - y(1);
u = u - u(1);

n = 1; % Output order
m = 2; % Input order (2 terms for B_d)
lambda = 0.99; % Forgetting factor (adjust between 0.98-1 for noise rejection)

N = length(y); % Number of samples
theta_hist = zeros(3, N);

% Initialize parameters
theta_hat = zeros(n + m, 1); % [a1, b0, b1]
P = 1e6 * eye(n + m);  % Large initial covariance

for k = 3:N
    x_k = [y(k-1); u(k); u(k-1)];  % Include past inputs
    
    y_hat = x_k' * theta_hat;
    
    K = (P * x_k) / (lambda + x_k' * P * x_k);  % Gain update with forgetting factor
    theta_hat = theta_hat + K * (y(k) - y_hat); % Update parameter estimate
    for i = 1:length(theta_hat)
        theta_hist(i,k) = theta_hat(i);
    end
    P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k); % Covariance update
end

% Display estimated parameters
disp('Final Estimated Parameters (theta_hat):');
disp(theta_hat);

% Assign estimated coefficients
a1 = theta_hat(1);  % y(k-1) coefficient
b0 = theta_hat(2);  % u(k) coefficient
b1 = theta_hat(3);  % u(k-1) coefficient

% Predicted output using estimated parameters
yh = zeros(N, 1);
for i = 2:N-1
    yh(i+1) = a1 * yh(i) + b0 * u(i) + b1 * u(i-1);
end

figure;
plot(y, 'r', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(yh, 'b--', 'DisplayName', 'Identificación', LineWidth=1.8);
legend;
title('Comparación de salida, método recursivo');
xlabel('Tiempo');
ylabel('Salida');
fontsize(16, "points");

figure;
plot(theta_hist(1,1:end), 'r', 'DisplayName', 'a1', LineWidth=1.8);
hold on;
plot(theta_hist(2,1:end), 'b--', 'DisplayName', 'b0', LineWidth=1.8);
hold on;
plot(theta_hist(3,1:end), 'g--', 'DisplayName', 'b1', LineWidth=1.8);
legend;
title('parameter update');
xlabel('Tiempo');
ylabel('Parameters');
fontsize(16, "points");

