% Method of least squares, direct

y = out.output_simout; % Output data
u = out.input_simout;  % Input data

% Subtract the initial output and input
y = y-y(1);
u = u-u(1);

n = 1; % Output order
m = 1; % Input order

N = length(y); % Number of samples

% Preallocate Phi and Y
Phi = zeros(N-1, n + m);
Y = zeros(N-1, 1);

% Construct Phi and Y matrices
for k = 2:N
    phi_y = y(k-1);
    phi_u = u(k-1);
    Phi(k-1, :) = [phi_y, phi_u];
    Y(k-1) = y(k);
end

theta_hat = (Phi' * Phi) \ (Phi' * Y);

disp('Estimated Parameters (theta_hat):');
disp(theta_hat); 
%Comprobación
a = theta_hat(1);
b = theta_hat(2);

yh=zeros(N,1);
for i=1:N-1
    yh(i+1) = a*yh(i)+b*u(i);
end

figure;
plot(y, 'r', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(yh, 'b--', 'DisplayName', 'Identificación', LineWidth=1.8);
legend;
fontsize(16, "points");
title("Comparación de salida, método directo");
xlabel('Tiempo');
ylabel('Salida');

