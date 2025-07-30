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

k_list = [0.05, 0.15, 0.30, 0.50, 0.90];

% ZOH
Ts = 1; % Tiempo del muestreo

G = tf([0 dfdk_eval], [1 -dfdh_eval]);
SC = ss(G);
Gz = c2d(G, Ts);
%Gz_dist = c2d(G, Ts*5);
%[B,A] = tfdata(Gz, 'v');
%[B_dist,A_dist] = tfdata(Gz_dist, 'v');

[stepG, tOutG] = step(G);
[stepGz, tOutGz] = step(Gz);

figure;
stairs(tOutGz, stepGz, 'b', 'DisplayName', 'Real', LineWidth=1.8)
hold on;
plot(tOutG, stepG, 'r', 'DisplayName', 'Real', LineWidth=1.8);
grid;
legend('G','Gz','Location','SouthWest');
xlabel("Time (seconds)");
fontsize(16, "points");
title('Step response');


n = 300;
plc = [.86 .86];
Alc = poly(plc);
alc1 = Alc(2);
alc2 = Alc(3);

prec = Ts/10;
t = 0:prec:Ts;
yc = zeros(size(t));
uu = [];
yy = [];
vv = [];
x = zeros(1,length(SC.a));
y = zeros(n,1);
u = zeros(n,1);
r = ones(n,1)*0.6; % Constant reference

%Simulación

vc = zeros(size(t))'; % Disturbance

% Inicializar RLS
theta_hat = [-2; 2; -0.2]; % Better initial guess for [s0; s1; t0]
P = 10 * eye(3); % 3x3 covariance matrix
lambda = 0.90; % Balanced forgetting factor

theta_hist = zeros(3, n);

% Reference model
plc = [0.86, 0.86];
Alc = poly(plc); % Alc = [1, -1.72, 0.7396]
y_m = zeros(n, 1);
for i = 3:n
    % Use precomputed r(i-1) (already includes the step change)
    if i > 200; r(i) = 0.2; end
    y_m(i) = -Alc(2)*y_m(i-1) - Alc(3)*y_m(i-2) + (1 - (-Alc(2)-Alc(3))) * r(i-1);
end

% Plot reference model
figure;
plot(y_m, 'r', 'DisplayName', 'y_m', LineWidth=1.8);
hold on;
plot(r, 'g--', 'DisplayName', 'Referencia', LineWidth=1.8);
legend('Location', 'southeast');
ylim([0 0.65]);
xlabel('Tiempo (s)');
fontsize(16, "points");

for i = 5:n
    y(i) = yc(end); % Output from previous control action

    % Regression vector with integral term
    x_k = [-y(i); -y(i-1); r(i)];
    
    % RLS update
    K = (P * x_k) / (lambda + x_k' * P * x_k);
    theta_hat = theta_hat + K * (y(i) - y_m(i));
    P = (P - K * x_k' * P) / lambda;
    
    % Extract parameters
    s0 = theta_hat(1);
    s1 = theta_hat(2);
    t0 = theta_hat(3);
    
    % Incremental control law with saturation
    u(i) = u(i-1) + t0 * r(i) - s0 * y(i) - s1 * y(i-1);
    %u_max = 1.0;
    %u_min = -0.5;
    %u(i) = max(min(u(i), u_max), u_min);
    
    uc = u(i)*ones(size(t))'; % Bloqueador actuacion

    if i > n/3 % Se introduca perturbacion
        vc = -.2*ones(size(t))';
    end
    
    % Simular el sistema
    [yc,ts,x] = lsim(SC,uc+vc,t,x(length(x),:));
    uu = [uu; uc(2:length(uc))];
    yy = [yy; yc(2:length(yc))];
    vv = [vv; vc(2:length(yc))];
end

num_samples = length(yy);
time = (0:0.1:(num_samples-1)/10);

figure;
plot(theta_hist(1,1:end), 'r', 'DisplayName', 's_0', LineWidth=1.8);
hold on;
plot(theta_hist(2,1:end), 'b', 'DisplayName', 's_1', LineWidth=1.8);
legend('Location', 'east');
ylabel('Valor')
xlabel('Tiempo (s)');
fontsize(16, "points");

figure;
plot(r, 'g--', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, uu, 'b', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, yy, 'r', 'DisplayName', 'Real', LineWidth=1.8)
grid;
legend('Referencia','Control input','Output','Location','SouthEast');
xlabel('Tiempo (s)');
fontsize(16, "points");