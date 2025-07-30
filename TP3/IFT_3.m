% Initial PID Parameters
Kp = 20; Ki = 1; Kd = 50;
parameters = [Kp; Ki; Kd];
gamma = 10; % Learning rate
lambda = 0.1; % Regularization term
num_iterations = 100; % Number of iterations
%delta = 0.1; % Perturbation size
%Kp_old = zeros(1, num_iterations);

for iter = 1:num_iterations
    fprintf('Iteration %d\n', iter);

    % Run Simulink model
    Kp = parameters(1);
    Ki = parameters(2);
    Kd = parameters(3);
    out = sim('practico3_IFT.slx');
    y = out.nivout.Data;
    u = out.ctrout.Data;
    t = out.nivout.Time;

    % Cost function
    y_d = h_list(4) * ones(size(t));
    e = y - y_d;
    J = 10*mean(e.^2) + lambda * mean(u.^2);

    % Gradient approximation
    dJ = zeros(1, 3);
    for j = 1:3
        parameters_p = parameters;
        delta = 0.01*abs(parameters_p(j));
        parameters_p(j) = parameters(j) + delta;
        Kp = parameters_p(1);
        Ki = parameters_p(2);
        Kd = parameters_p(3);
        %fprintf("Kp = %.4f Ki = %.4f Kd = %.4f \n", parameters_p(1), parameters_p(2), parameters_p(3));
        out_perturbed = sim("practico3_IFT.slx");
        % Perturbed cost function
        y_perturbed = out_perturbed.nivout.Data;
        u_perturbed = out_perturbed.ctrout.Data;
        e_perturbed = y_perturbed - y_d;
        J_perturbed = 10*mean(e_perturbed.^2) + lambda * mean(u_perturbed.^2);
        dJ(j) = (J_perturbed-J) / delta;
    end
    %Kp = Kp + delta;
    
    out_perturbed = sim('practico3_IFT.slx');

    % gradient
    %dJ = (J_perturbed - J) / delta;
    % Update parameters with constraints
    parameters = parameters - gamma*dJ;
    % Current PID parameters
    
    % Log progress
    fprintf('Iteration %d: J = %.4f, Kp = %.4f, Ki = %.4f, Kd = %.4f\n', ...
            iter, J, parameters(1), parameters(2), parameters(3));
end

% Display final parameters
disp('Final PID Parameters:');