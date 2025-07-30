clear all;

r_0 = 0.05; %Radio menor
r_o = 0.0254; %Radio de orifício de salida
g = 9.81;
theta = 60;

a = pi*r_o^2; %La sección de tubería de salida

F_in = 0.0035;

%k_out = 1; %Porcentaje de apertura

h_bar = 0.2; %punto de trabajo
k_bar = F_in / (a*sqrt(2*9.81*h_bar)); %k en el punto de trabajo

F_out = a*sqrt(2*g);

%%Linearización
syms h k

F_out_lin = k*a*sqrt(2*g);

f = (F_in - F_out_lin*sqrt(h)) / (pi * ((h^2)/3 + 2*r_0*h/sqrt(3) + r_0^2));

dfdh = diff(f, h);
dfdh_eval = double(vpa(subs(dfdh, [h, k], [h_bar, k_bar])));

dfdk = diff(f, k);
dfdk_eval = double(vpa(subs(dfdk, [h, k], [h_bar, k_bar])));

% Discrete system
Ts = 1;  % Sampling time
G = tf([0 dfdk_eval], [1 -dfdh_eval], 'InputDelay', 0);
Gz = c2d(G, Ts, 'zoh');
[B,A] = tfdata(Gz, 'v');

% STR parameters
n = 300;            % Simulation steps
lambda = 0.95;      % Forgetting factor
na = 2; nb = 2;     % System order
np = na + nb;       % Total parameters

% Pole placement parameters
plc = [0.3 0.3 0.3]; % Desired closed-loop poles
Alc = poly(plc);      % Characteristic polynomial

% Initialization
Aest = 0.5*ones(n, np);   % Parameter estimates (more robust init)
P = 1000*eye(np);         % Covariance matrix
th = [0.5; 0.5; 0.5; 0.5];% Initial parameters (non-zero)
xi = zeros(np,1);         % Regressor vector

% Integral action parameters
x1 = -0.7;               % Integral filter coefficient
X = [1 x1];              % Integral polynomial

% Simulation buffers
uu = zeros(n,1);
yy = zeros(n,1);
r = 0.6*ones(n,1);       % Reference signal
vc = 0;                  % Disturbance

% Main loop
for i = 5:n
    % Update reference
    if i > 2*n/3
        r(i) = 0.2;      % Reference change
    end
    
    % System output
    yy(i) = yy(i-1) + B(2)*uu(i-1) - A(2)*yy(i-1);  % Simulated output
    
    % Regressor vector
    xi = [-yy(i-1); -yy(i-2); uu(i-1); uu(i-2)];
    
    % RLS estimation
    y_hat = xi'*th;
    eps = yy(i) - y_hat;
    K = P*xi/(lambda + xi'*P*xi);
    th = th + K*eps;
    P = (P - K*xi'*P)/lambda;
    Aest(i,:) = th';
    
    % Extract parameters with projection
    th = max(min(th, 2), 0.1);  % Keep parameters bounded
    
    % Pole placement with Sylvester matrix
    A_est = [1 th(1) th(2)];  % Degree 2 (1 + 2 coefficients)
    B_est = [0 th(3) th(4)];  % Degree 2 (0 + 2 coefficients)
    
    % Build proper Sylvester matrix (4x4 for 2nd order system)
    MA = toeplitz([A_est zeros(1, length(A_est)-1)], ...
                 zeros(length(A_est)+length(B_est)-2, 1));
    MB = toeplitz([B_est zeros(1, length(B_est)-1)], ...
                 zeros(length(A_est)+length(B_est)-2, 1));
    M = [MA; MB];
    
    % Adjust desired polynomial to match dimensions
    desired_order = length(A_est) + length(B_est) - 1;
    Alc_adj = Alc(1:desired_order);
    
    % Solve for controller polynomials
    RS = M\Alc_adj';  % Now M is 4x4 and Alc_adj is 1x4
    Rr = RS(1:length(A_est)-1)';
    Sr = RS(length(A_est):end)';
    
    % Add integral action
    aux1 = X*sum(B_est);
    aux2 = conv(aux1, Rr);
    aux3 = conv([1 x1], sum(Rr));
    aux4 = conv(aux3, B_est);
    R = aux2 - [0 aux4];
    S = conv(aux1, Sr) + conv(aux3, A_est);
    
    % Control law
    uu(i) = (sum(R)*r(i) - S(1)*yy(i) - S(2)*yy(i-1))/sum(R);
    
    % Anti-windup and saturation
    uu(i) = max(min(uu(i), 1), -1);  % ±1 m³/s limits
end

% Plot results
figure;
subplot(2,1,1);
plot(yy, 'b', 'LineWidth',1.5); hold on;
plot(r, 'k--', 'LineWidth',1.5);
title('System Output');
legend('Actual','Reference');

subplot(2,1,2);
plot(Aest(:,1:2), 'LineWidth',1.5); hold on;
plot(Aest(:,3:4), '--', 'LineWidth',1.5);
title('Parameter Estimates');
legend('a1','a2','b1','b2');