% Recursive Least Squares (RLS) - Adjusted for STR

y = out.nivout; % Output data
u = out.ctrout; % Input data

% Subtract the initial output and input
y = y - y(1);
u = u - u(1);

n = 1; % Output order
m = 2; % Input order (we now estimate two input terms)

N = length(y); % Number of samples

% Initialize parameters
theta_hat = zeros(n + m, 1); % Now has 3 parameters
P = 1e6 * eye(n + m);        % Large initial covariance matrix

for k = 3:N  % Start from 3 because we need u(k-1) and u(k-2)
    x_k = [y(k-1); u(k); u(k-1)];  % Include past input
    
    y_hat = x_k' * theta_hat;
    
    theta_hat = theta_hat - P * x_k * (y_hat - y(k));
    
    P = P - (P * (x_k * x_k') * P) / (1 + x_k' * P * x_k);
    
end

% Display final estimated parameters
disp('Final Estimated Parameters (theta_hat):');
disp(theta_hat);

% Assign values
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