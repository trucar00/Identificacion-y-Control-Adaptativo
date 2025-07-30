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
[B,A] = tfdata(Gz, 'v');
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
plc = .86;
Alc = poly([plc plc]);
alc1 = Alc(2);
alc2 = Alc(3);

a1 = A(2);
b1 = B(2);

s0 = 1/b1*(alc1-a1+1);
s1 = 1/b1*(alc2+a1);
t0 = sum(Alc)/b1;



prec = Ts/10;
t = 0:prec:Ts;
yc = zeros(size(t));
uu = [];
yy = [];
vv = [];
%yym = [];
%tt = [];
x = zeros(1,length(SC.a));
%lx = length(SC.a);
y = zeros(n,1);
%yd = zeros(n,1);
u = zeros(n,1);
%u_ref = zeros(n,1);
%u_dist = zeros(n,1);
r = ones(n,1)*0.6; % Constant reference
%rf = zeros(n,1);

%ym = zeros(n,1);
%Simulación

vc = zeros(size(t))'; % Disturbance

for i = 5:n
    if i>n*2/3; r(i) = 0.2; end
    y(i) = yc(end); % Sample of the output

    %u_ref(i) = t0*r(i) - s0*y(i); % Seguir la referencia
    %u_dist(i) = -(vc(end)); % Corregir perturbaciones

    %uc = (u_ref(i)+u_dist(i))*ones(size(t))'; % Bloqueador actuacion

    %cálculo de la accion de control
    u(i) = u(i-1)+t0*r(i)-s0*y(i)-s1*y(i-1);

    % bloqueador actuación
    uc = u(i)*ones(size(t))';

    if i > n/3 % Introduce disturbance
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
plot(r, 'g--', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, uu, 'b', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, yy, 'r', 'DisplayName', 'Real', LineWidth=1.8)
grid;
legend('Referencia','Control input','Output','Location','East');
xlabel('Tiempo (s)');
fontsize(16, "points");

%time = (0:n-1);
%u_timeseries = timeseries(u, time);