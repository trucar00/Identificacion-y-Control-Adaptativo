% Construir la theta con los parámetros obtenidos desde Ziegler-Nichols
Kp = 8.8800; Ki = 0.2916; Kd = 67.5990;
theta_IFT = [Kp, Ki, Kd];

gamma = 1; % Tasa de aprendizaje
alpha = 100; % Penalización, referencia
lambda = 0.1; % Penalización, control

num_iteraciones = 150;

for iter = 1:num_iteraciones

    % Simular el modelo
    out = sim('practico3_IFT.slx');
    y = out.nivout.Data;
    u = out.ctrout.Data;

    % Función de costos
    y_d = h_list(3) * ones(size(y));
    J = (100*sum((y - y_d).^2) + lambda*sum(u.^2))/(2*length(y));

    % Aproximacion del gradiente
    gradient = zeros(1, 3);
    for j = 1:3
        theta_perturbado = theta_IFT;
        delta = 0.1 * abs(theta_IFT(j)); % La delta depende del tamaño del parámetro
        theta_perturbado(j) = theta_perturbado(j) + delta;

        % Simular el modelo con parámetros perturbado
        Kp = theta_perturbado(1);
        Ki = theta_perturbado(2);
        Kd = theta_perturbado(3);
        out_perturbado = sim('practico3_IFT.slx');

        % Funcion de costos perturbado
        y_perturbado = out_perturbado.nivout.Data;
        u_perturbado = out_perturbado.ctrout.Data;
        J_perturbado = (100*sum((y_perturbado-y_d).^2) + lambda*sum(u_perturbado.^2))/(2*length(y_perturbado));

        % Calcular el gradiente
        gradient(j) = (J_perturbado - J) / delta;
    end

    % Actualizar los parámetros
    theta_IFT = theta_IFT - gamma * gradient;
    Kp = theta_IFT(1);
    Ki = theta_IFT(2);
    Kd = theta_IFT(3);

    % Log progress
    fprintf('Iteración %d: J = %.4f, Kp = %.4f, Ki = %.4f, Kd = %.4f\n', ...
            iter, J, theta_IFT(1), theta_IFT(2), theta_IFT(3));
end

% Display final parameters
disp('Final PID Parameters:');
fprintf('Kp = %.4f, Ki = %.4f, Kd = %.4f\n', theta_IFT(1), theta_IFT(2), theta_IFT(3));