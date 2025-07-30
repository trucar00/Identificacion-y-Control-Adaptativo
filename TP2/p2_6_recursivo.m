% Recursive Least Squares (RLS)

y = out.nivout; % Output data
u = out.ctrout;  % Input data

% Subtract the initial output and input
y = y - y(1);
u = u - u(1);

n = 1; % Output order
m = 1; % Input order

N = length(y); % Number of samples

% Initialize parameters
theta_hat = zeros(n + m, 1); % Initial parameter estimations
P = 1e6 * eye(n + m);        % Large initial covariance matrix

for k = 2:N
    x_k = [y(k-1); u(k-1)];
    
    y_hat = x_k' * theta_hat;
    
    theta_hat = theta_hat - P * x_k * (y_hat - y(k));
    
    P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k);
    
end

% Display final estimated parameters
disp('Final Estimated Parameters (theta_hat):');
disp(theta_hat);

% Compare estimated output with actual output
a = theta_hat(1);
b = theta_hat(2);

yh = zeros(N, 1); % Predicted output
for i = 1:N-1
    yh(i+1) = a * yh(i) + b * u(i);
end

% Plot results
figure;
plot(y, 'r', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(yh, 'b--', 'DisplayName', 'Identificación', LineWidth=1.8);
legend;
title('Comparación de salida, método recursivo');
xlabel('Tiempo');
ylabel('Salida');
fontsize(16, "points");