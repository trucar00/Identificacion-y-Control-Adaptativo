clear all;

r_0 = 0.05; %Radio menor
r_o = 0.0254; %Radio de orifício de salida
g = 9.81;
theta = 60;

a = pi*r_o^2; %La sección de tubería de salida

F_in = 0.0035;

%k_out = 1; % Porcentaje de apertura

h_bar = 0.2; % Punto de trabajo
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

%%Discretización
Ts = 1; %tiempo del muestreo

A_d = exp(dfdh_eval*Ts);
B_d = dfdh_eval^(-1) * (A_d-1) * dfdk_eval;


% ZOH
G = tf([0 dfdk_eval], [1 -dfdh_eval]);
SC = ss(G);
Gz = c2d(G, Ts);
[Bd,Ad] = tfdata(Gz, 'v');
a1 = Ad(2);
b1 = Bd(2);

[stepG, tOutG] = step(G);
[stepGz, tOutGz] = step(Gz);

figure;
stairs(tOutGz, stepGz, 'b', 'DisplayName', 'Real', LineWidth=1.8)
hold on;
plot(tOutG, stepG, 'r', 'DisplayName', 'Real', LineWidth=1.8);
grid;
legend('G(z)','G','Location','NorthEast');
xlabel("Time (s)");
fontsize(16, "points");
%title('Step response');

n = 300;
plc = .86; % Polo deseado
Alc = poly(plc);
alpha = Alc(2);
s0 = (alpha-a1) / b1;
t0 = (1+alpha) / b1;
prec = Ts/10;
t = 0:prec:Ts;
yc = zeros(size(t));
uu = [];
yy = [];
vv = [];
x = zeros(1,length(SC.a));
y = zeros(n,1);
u = zeros(n,1);
r = ones(n,1)*0.6; % Referencia constante

vc = zeros(size(t))'; % Perturbacion

for i = 5:n
    if i>n*2/3; r(i) = 0.2; end
    y(i) = yc(end); % Muestreo de la salida

    u(i) = t0*r(i)-s0*y(i); % Accion de control
    
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
plot(r, 'g--', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, uu, 'b', 'DisplayName', 'Real', LineWidth=1.8);
hold on;
plot(time, yy, 'r', 'DisplayName', 'Real', LineWidth=1.8)
grid;
legend('Referencia','Acción de control','Salida','Location','East');
xlabel('Tiempo (s)');
fontsize(16, "points");

%time = (0:n-1);
%u_timeseries = timeseries(u, time);