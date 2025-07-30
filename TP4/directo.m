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
Ts = 10; % Tiempo del muestreo

G = tf([0 dfdk_eval], [1 -dfdh_eval]);
SC = ss(G);
Gz = c2d(G, Ts);
%Gz_dist = c2d(G, Ts*5);
%[B,A] = tfdata(Gz, 'v');
%[B_dist,A_dist] = tfdata(Gz_dist, 'v');


n = 300;
plc = [.6 .6];
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
P = 10*eye(2);             
theta_hat = [-0.3; -3]; % Se puede "tune"
lambda = 0.85; % Factor de olvido

theta_hist = zeros(2, n);

alpha = 0.60; % Polo deseado
y_m = zeros(n, 1); % Modelo

for i = 5:n
    if i > 200; r(i) = 0.2; end
    y_m(i) = alpha * y_m(i-1) + (1 - alpha) * r(i-1);

    y(i) = yc(end); % Muestreo de la salida

    % Metodo de minimos cuadrados
    x_k = [-y(i-1); r(i-1)];
    K = (P * x_k) / (lambda + x_k' * P * x_k);
    theta_hat = theta_hat + K * (y(i) - y_m(i));
    P = (P - K * x_k' * P) / lambda;
    
    theta_hist(:, i) = theta_hat;
    
    % Parametros de la ley de control
    s0 = theta_hat(1);
    t0 = theta_hat(2);
    
    u(i) = u(i-1) + t0 * r(i) - s0 * y(i);
    % Can I change this to: u(i) = u(i-1) + t0 * r(i) - s0 * y(i); to get
    % incremental control?
    
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
plot(y_m, 'r', 'DisplayName', 'y_m', LineWidth=1.8)
hold on;
plot(r, 'g--', 'DisplayName', 'Referencia', LineWidth=1.8);
legend('Location', 'southeast');
ylim([0 0.65])
xlabel('Tiempo (s)');
fontsize(16, "points");

figure;
plot(theta_hist(1,1:end), 'r', 'DisplayName', 's_0', LineWidth=1.8);
hold on;
plot(theta_hist(2,1:end), 'b', 'DisplayName', 't_0', LineWidth=1.8);
legend('Location', 'southeast');
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