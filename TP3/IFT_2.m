% Initialize Parameters
K = [1; 1; 1]; % Initial gains for Kp, Ki, Kd
gamma = 0.1; % Gradient descent step size
lambda = 0.1; % Weighting for control cost
beta = ones(3, 1); % Weighting for partial derivatives
perturbation = 1e-3; % Perturbation size (e.g., 0.1% of the gain)

% Desired output (reference signal)
yd = h_list(4) * ones(20001, 1); % Example reference signal (size matches y)

% Iterative Optimization
for iter = 1:100
    % Simulate the system for nominal gains
    Kp = K(1);
    Ki = K(2);
    Kd = K(3);
    out_nominal = sim("practico3_IFT.slx");

    % Extract nominal outputs
    y_nominal = out_nominal.nivout.Data;
    u_nominal = out_nominal.ctrout.Data;

    % Preallocate perturbed outputs
    y_perturbed = zeros(length(y_nominal), 3); % Columns for Kp, Ki, Kd
    u_perturbed = zeros(length(u_nominal), 3);

    % Perturb each gain and simulate
    for gain_idx = 1:3
        % Perturb the current gain
        K_perturbed = K;
        K_perturbed(gain_idx) = K_perturbed(gain_idx) + perturbation;

        % Update Simulink parameters
        Kp = K_perturbed(1);
        Ki = K_perturbed(2);
        Kd = K_perturbed(3);

        % Run the simulation with perturbed gains
        out_perturbed = sim("practico3_IFT.slx");

        % Extract perturbed outputs
        y_perturbed(:, gain_idx) = out_perturbed.nivout.Data;
        u_perturbed(:, gain_idx) = out_perturbed.ctrout.Data;
    end

    % Combine nominal and perturbed signals
    y = [y_nominal, y_perturbed];
    u = [u_nominal, u_perturbed];

    % Update the PID gains using the custom gradient function
    K = pidift(K, yd, y, u, gamma, lambda, beta);

    % Log or plot the results (optional)
    fprintf('Iteration %d: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', iter, K(1), K(2), K(3));

    
    %K_old = K;
end

function K = pidift(K, yd, y, u, gamma, lambda, beta)

    K = K(:);
    
    if nargin < 4
        gamma = 0.1;
    end
    if nargin < 5
        lambda = 0;
    end
    if nargin < 6
        beta = ones(3,1);
    end
    
    % Extract Kp and Ki, Kd is not used directly in the cost calculations.
    Kp = K(1);
    Ki = K(2);
    
    % Extract the test signals.
    y1 = y(:,1);
    y2 = y(:,2);
    y3 = y(:,3);
    u1 = u(:,1);
    u2 = u(:,2);
    u3 = u(:,3);
    
    % Construct the gradient PID controller functions. These are the partial
    % derivatives of the PID controller TF.
    dKp = tf([1, -1], [Kp + Ki, -Kp], 1);
    dKi = tf([1, 0], [Kp + Ki, -Kp], 1);
    dKd = tf([1, -2, 1], [Kp + Ki, -Kp, 0], 1);
    
    % Calculate the gradient signals.
    dydKp = beta(1)*lsim(dKp, y2);
    dydKi = beta(2)*lsim(dKi, y2);
    dydKd = beta(3)*lsim(dKd, y2 - y3);
    
    dudKp = lsim(dKp, u2);
    dudKi = lsim(dKi, u2);
    dudKd = lsim(dKd, u2 - u3);
    
    % Calculate the gradient of the cost function for each gain.
    dJ = 0;
    for i = 1:length(y1)
        dJ = dJ + (y1(i) - yd(i)) * [dydKp(i); dydKi(i); dydKd(i)] ...
            + lambda * u1(i) * [dudKp(i); dudKi(i); dudKd(i)];
    end
    dJ = dJ/length(y1);
    
    % Calculate the Guass-Newton-ish gradient matrix.
    R = 0;
    for i = 1:length(y1)
        dydrho = [dydKp(i); dydKi(i); dydKd(i)];
        dudrho = [dudKp(i); dudKi(i); dudKd(i)];
        R = R + dydrho*dydrho' + dudrho*dudrho';
    end
    R = R/length(y1);
    
    % Calculate the new gains.
    K = K - gamma*inv(R)*dJ;
end
